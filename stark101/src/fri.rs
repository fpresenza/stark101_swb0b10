use lambdaworks_math::traits::{AsBytes, ByteConversion};
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
};
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_crypto::merkle_tree::{
    merkle::MerkleTree,
    backends::types::Keccak256Backend, 
    proof::Proof
};
use lambdaworks_crypto::fiat_shamir::{
    is_transcript::IsTranscript,
    default_transcript::DefaultTranscript
};
use lambdaworks_math::unsigned_integer::element::U256;

use crate::poly;

pub struct InclusionProof<F: IsField> (
    pub FieldElement<F>,
    pub Proof<[u8; 32]>
);

#[derive(Clone)]
pub struct ValidationData<F: IsField> {
    proof: Proof<[u8; 32]>,
    sym_eval: FieldElement<F>,
    sym_proof: Proof<[u8; 32]>,
}

#[derive(Clone)]
pub struct FriLayer<F: IsField> {
    root: [u8; 32],
    queries: Vec<ValidationData<F>>,
}

pub fn trace_inclusion_proofs<F>(
        query_indices: &Vec<usize>,
        domain_size: usize,
        trace_poly_eval: &[FieldElement<F>],
        trace_poly_tree: &MerkleTree<Keccak256Backend<F>>,
        transcript: &mut DefaultTranscript<F>,
    ) -> Vec<[InclusionProof<F>; 3]> 
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    query_indices
        .iter()
        .map(|i|{
            let idx = i.to_owned();
            let idx1 = (idx + 1) % domain_size;
            let idx2 = (idx + 2) % domain_size;
            transcript.append_bytes(&idx.to_be_bytes());
            
            [
            InclusionProof(trace_poly_eval[idx].to_owned(), trace_poly_tree.get_proof_by_pos(idx).unwrap()),
            InclusionProof(trace_poly_eval[idx1].to_owned(), trace_poly_tree.get_proof_by_pos(idx1).unwrap()),
            InclusionProof(trace_poly_eval[idx2].to_owned(), trace_poly_tree.get_proof_by_pos(idx2).unwrap())
            ]
        })
        .collect()
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        query_indices: Vec<usize>,
        transcript: &mut DefaultTranscript<F>
    ) -> Vec<FriLayer<F>>
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let mut polynomial = polynomial.clone();
    let mut offset = offset.clone();
    let number_of_foldings = (usize::BITS - polynomial.degree().leading_zeros()) as usize;
    let mut fri_layers = Vec::<FriLayer<F>>::with_capacity(number_of_foldings + 1);

    // commit to evaluations
    let (eval, tree) = commit(&polynomial, domain_size, &offset);
    transcript.append_bytes(&tree.root);
    println!(
        "Layer 0: \n \t Appending root of composition polynomial (degree {:?}) to transcript.",
         polynomial.degree()
    );

    // Generate inclusion proofs, validation data and append to layer
    fri_layers.push(
        FriLayer {
            root: tree.root,
            queries: query_indices.iter().map(|i| { 
                let idx = i.to_owned();
                let sym_idx = (idx + domain_size / 2) % domain_size;
        
                ValidationData {
                    proof: tree.get_proof_by_pos(idx).unwrap(),
                    sym_eval: eval[sym_idx].to_owned(),
                    sym_proof: tree.get_proof_by_pos(sym_idx).unwrap()
                }
            })
            .collect::<Vec<ValidationData<F>>>()
        }
    );

    // recursive foldings
    for layer in 1..=number_of_foldings {
        let beta = transcript.sample_field_element();
        (polynomial, domain_size, offset) = fold(polynomial, domain_size, offset, beta);

        let (eval, tree) = commit(&polynomial, domain_size, &offset);
        transcript.append_bytes(&tree.root);
        println!(
            "Layer {:?}: \n \t Appending root of folded polynomial (degree {:?}) to transcript.",
            layer,
            polynomial.degree()
        );

        // append layer
        fri_layers.push(
            FriLayer {
                root: tree.root,
                queries: query_indices.iter().map(|i| { 
                    let idx = i.to_owned() % domain_size;
                    let sym_idx = (idx + domain_size / 2) % domain_size;
        
                    ValidationData {
                        proof: tree.get_proof_by_pos(idx).unwrap(),
                        sym_eval: eval[sym_idx].to_owned(),
                        sym_proof: tree.get_proof_by_pos(sym_idx).unwrap()
                    }
                })
                .collect::<Vec<ValidationData<F>>>()
            }
        );
    }

    fri_layers
}

fn commit<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> (Vec<FieldElement<F>>, MerkleTree<Keccak256Backend<F>>)
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let eval = Polynomial::evaluate_offset_fft::<F>(
        polynomial,
        1, 
        Some(domain_size),
        offset
    ).unwrap();
    
    let tree = MerkleTree::<Keccak256Backend<F>>::build(&eval);

    (eval, tree)
}

fn fold<F: IsField>(
        polynomial: Polynomial<FieldElement<F>>,
        domain_size: usize,
        offset: FieldElement<F>,
        beta: FieldElement<F>
    ) -> (Polynomial<FieldElement<F>>, usize, FieldElement<F>) {
    (poly::fold_polynomial(&polynomial, &beta),
    domain_size / 2,
    offset.square())
}

pub fn sample_queries<F>(
        num_queries: usize,
        domain_size: usize,
        transcript: &mut DefaultTranscript<F>
    ) -> Vec<usize> 
    where 
        F: IsField,
        FieldElement<F>: AsBytes + ByteConversion {

        (0..num_queries)
        .map(|_| {
            let query_index = U256::from_bytes_be(&transcript.sample()).unwrap();
            let(_, query_index) = query_index.div_rem(&U256::from(domain_size as u64));
            query_index.limbs[3] as usize
        })
        .collect::<Vec<usize>>()
}