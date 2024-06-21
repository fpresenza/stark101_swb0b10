use lambdaworks_math::unsigned_integer::element::U256;
use lambdaworks_math::field::{
	traits::IsField,
    element::FieldElement
};

use crate::fri::FriLayer;

#[derive(Clone)]
pub struct PublicInput<F: IsField> (
	pub U256,
	pub usize,
	pub usize,
	pub usize,
	pub FieldElement<F>,
	pub FieldElement<F>
);

#[derive(Clone)]
pub struct StarkProof<F: IsField> {
	pub trace_root: [u8; 32],
	pub fri_layers: Vec<FriLayer<F>>
}

