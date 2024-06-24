use lambdaworks_math::field::traits::IsPrimeField;
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

use crate::poly;

pub type FriCommitment<F> = Vec<FriLayer<F>>;

#[derive(Clone)]
pub struct ValidationData<F: IsField> {
    pub proof: Proof<[u8; 32]>,
    pub sym_eval: FieldElement<F>,
    pub sym_proof: Proof<[u8; 32]>,
}

#[derive(Clone)]
pub struct FriLayer<F: IsField> {
    pub root: [u8; 32],
    pub validation_data: Vec<ValidationData<F>>,
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        query_indices: Vec<usize>,
        transcript: &mut DefaultTranscript<F>
    ) -> Vec<FriLayer<F>>
    where
        F: IsField + IsFFTField + IsPrimeField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let mut polynomial = polynomial.clone();
    let mut offset = offset.clone();
    let number_of_foldings = (usize::BITS - polynomial.degree().leading_zeros()) as usize;
    let mut fri_layers = Vec::<FriLayer<F>>::with_capacity(number_of_foldings + 1);

    // commit to evaluations
    let (eval, tree) = commit(&polynomial, domain_size, &offset);
    transcript.append_bytes(&tree.root);

    // Generate inclusion proofs, validation data and append to layer
    fri_layers.push(
        FriLayer {
            root: tree.root,
            validation_data: query_indices.iter().map(|i| { 
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
    for _ in 1..=number_of_foldings {
        let beta = transcript.sample_field_element();

        (polynomial, domain_size, offset) = fold(polynomial, domain_size, offset, beta);

        let (eval, tree) = commit(&polynomial, domain_size, &offset);
        transcript.append_bytes(&tree.root);

        // append layer
        fri_layers.push(
            FriLayer {
                root: tree.root,
                validation_data: query_indices.iter().map(|i| { 
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

pub fn decommit_and_fold<F>(
        layers: &[FriLayer<F>],
        domain_size: &usize,
        query_indices: &[usize],
        queries: &[FieldElement<F>],
        query_evals: &[FieldElement<F>],
        transcript: &mut DefaultTranscript<F>
    ) -> bool
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let mut domain_size = domain_size.to_owned();
    let mut queries = queries.to_owned();
    let mut query_evals = query_evals.to_owned();
    let mut sym_evals = Vec::<FieldElement<F>>::with_capacity(query_evals.len());

    // commit to evaluations
    let FriLayer{root, validation_data} = &layers[0];
    transcript.append_bytes(root);

    // verify first layer inclusion proofs and get next layer queries
    let num_queries = query_indices.len();
    for i in 0..num_queries {
        let idx = query_indices[i];
        let sym_idx = (idx + domain_size / 2) % domain_size;
        let eval = &query_evals[i];
        let ValidationData{proof, sym_eval, sym_proof} = &validation_data[i];
        sym_evals.push(sym_eval.clone());

        if !proof.verify::<Keccak256Backend<F>>(root, idx, eval) || !sym_proof.verify::<Keccak256Backend<F>>(root, sym_idx, sym_eval) {
            return false            
        }
    };

    // recursive foldings
    for layer in layers.iter().skip(1) {
        let beta = transcript.sample_field_element();
        
        domain_size /= 2;
        
        let FriLayer{root, validation_data} = layer;
        transcript.append_bytes(root);

        for i in 0..num_queries {
            query_evals[i] = curr_layer_query_evals(&queries[i], &query_evals[i], &sym_evals[i], &beta);
            queries[i] = queries[i].square();

            let idx = query_indices[i] % domain_size;
            let sym_idx = (idx + domain_size / 2) % domain_size;
            let eval = &query_evals[i];
            let ValidationData{proof, sym_eval, sym_proof} = &validation_data[i];
            sym_evals[i] = sym_eval.clone();

            if !proof.verify::<Keccak256Backend<F>>(root, idx, eval) || !sym_proof.verify::<Keccak256Backend<F>>(root, sym_idx, sym_eval) {
                return false            
            }
        }
    };

    true
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

pub fn curr_layer_query_evals<F: IsField>(
        query: &FieldElement<F>,
        eval: &FieldElement<F>,
        sym_eval: &FieldElement<F>,
        beta: &FieldElement<F>,
    ) -> FieldElement<F> {
    let query_inv = query.inv().unwrap();
    let two_inv = FieldElement::<F>::from(2_u64).inv().unwrap();
    ((eval + sym_eval) + beta * (eval - sym_eval) * query_inv) * two_inv
}