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

#[derive(Clone)]
pub struct Query<F: IsField> {
    index: Option<usize>,
    eval: Option<FieldElement<F>>,
    proof: Option<Proof<[u8; 32]>>,
    sym_eval: Option<FieldElement<F>>,
    sym_proof: Option<Proof<[u8; 32]>>,
}

#[derive(Clone)]
pub struct FriLayer<F: IsField> {
    root: [u8; 32],
    queries: Vec<Query<F>>,
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        num_queries: usize,
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

     // sample queries
    let query_indices = sample_queries(num_queries, domain_size, transcript);
    println!("\t Query indices: {:?}", query_indices);

    // get queries evaluations,add to transcript and generate inclusion proofs
    println!("\t Appending query indices to transcript.");
    // append layer
    fri_layers.push(
        FriLayer {
            root: tree.root,
            queries: query_indices.iter().map(|i| { 
                let idx = i.to_owned();
                let sym_idx = (idx + domain_size / 2) % domain_size;
                transcript.append_bytes(&idx.to_be_bytes());
        
                Query {
                    index: Some(idx),
                    eval: Some(eval[idx].to_owned()),
                    proof: Some(tree.get_proof_by_pos(idx).unwrap()),
                    sym_eval: Some(eval[sym_idx].to_owned()),
                    sym_proof: Some(tree.get_proof_by_pos(sym_idx).unwrap())
                }
            })
            .collect::<Vec<Query<F>>>()
        }
    );

    // recursive foldings
    for layer in 1..number_of_foldings {
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
        
                    Query {
                        index: None,
                        eval: None,
                        proof: Some(tree.get_proof_by_pos(idx).unwrap()),
                        sym_eval: Some(eval[sym_idx].to_owned()),
                        sym_proof: Some(tree.get_proof_by_pos(sym_idx).unwrap())
                    }
                })
                .collect::<Vec<Query<F>>>()
            }
        );
    }
    // final layer
    let beta = transcript.sample_field_element();
    (polynomial, domain_size, offset) = fold(polynomial, domain_size, offset, beta);

    let (eval, tree) = commit(&polynomial, domain_size, &offset);
    transcript.append_bytes(&tree.root);
    println!(
        "Layer {:?}: \n \t Appending root of folded polynomial (degree {:?}) to transcript.",
        number_of_foldings,
        polynomial.degree()
    );

    let constant_poly = polynomial.coefficients.first().unwrap();
    transcript.append_bytes(&constant_poly.to_bytes_be());
    println!("\t Appending constant polynomial to transcript.");

    // append layer
    fri_layers.push(
        FriLayer {
            root: tree.root,
            queries: query_indices.iter().map(|i| { 
                let idx = i.to_owned() % domain_size;
                let sym_idx = (idx + domain_size / 2) % domain_size;
        
                Query {
                    index: None,
                    eval: Some(constant_poly.to_owned()),
                    proof: Some(tree.get_proof_by_pos(idx).unwrap()),
                    sym_eval: Some(eval[sym_idx].to_owned()),
                    sym_proof: Some(tree.get_proof_by_pos(sym_idx).unwrap())
                }
            })
            .collect::<Vec<Query<F>>>()
        }
    );

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
    (
        poly::fold_polynomial(&polynomial, &beta),
        domain_size / 2,
        offset.square()
    )
}

fn sample_queries<F>(
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