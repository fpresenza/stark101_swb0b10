use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_crypto::merkle_tree::{
    merkle::MerkleTree,
    backends::types::Sha3_256Backend
};

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FE = FieldElement<F>;

// interpolation domain of size 1024 = 2^10
const INT_DOM_SIZE: usize = 0b10000000000;

// evaluation domain of size 8192 = 2^13 (blow-up factor is 2^3)
const EVAL_DOM_SIZE: usize = 0b10000000000000;
const BLOWUP_FACTOR: usize = 0b1000;

#[derive(Debug)]
pub struct Prover {}

impl Prover {

    pub fn generate_proof() {
        // ===================================
        // ==========|    Part 1:   |=========
        // === Statement, LDE & Commitment ===
        // ===================================
        // create vec to hold fibonacci square sequence
        let mut fib_sq = Vec::<FE>::with_capacity(INT_DOM_SIZE);
        fib_sq.push(FE::from(1_u64));
        fib_sq.push(FE::from(3141592_u64));
    
        for i in 2..INT_DOM_SIZE {
            let x = fib_sq[i-2];
            let y = fib_sq[i-1];
            fib_sq.push(x.pow(2_u64) + y.pow(2_u64));
        }
    
        // fft-interpolate the fibonacci square sequence
        let trace_poly = match Polynomial::interpolate_fft::<F>(&fib_sq) {
            Ok(poly) => poly,
            Err(e) => panic!("{:?}", e),
        };
        let int_dom_gen = F::get_primitive_root_of_unity(10_u64).unwrap();
        assert_eq!(trace_poly.coefficients.len(), INT_DOM_SIZE);
        assert_eq!(trace_poly.evaluate(&FE::from(1_u64)), FE::from(1_u64));
        assert_eq!(trace_poly.evaluate(&int_dom_gen), FE::from(3141592_u64));
        assert_eq!(
            trace_poly.evaluate(&int_dom_gen.pow(1022_u64)),
            FE::from_hex_unchecked("6A317721EF632FF24FB815C9BBD4D4582BC7E21A43CFBDD89A8B8F0BDA68252")
        );

        // fft-evaluate the fibonacci square sequence over a larger domain
        // of size (blow-up factor) * (interpolation domain size)
        // the offset is obtained as an outside not in the interpolation domain
        let offset = FE::from(2_u64);
        assert!(offset.pow(INT_DOM_SIZE as u64) != FE::from(1_u64));
        let trace_eval = match Polynomial::evaluate_offset_fft::<F>(
            &trace_poly, BLOWUP_FACTOR, None, &offset
        ) {
            Ok(poly) => poly,
            Err(e) => panic!("{:?}", e),
        };
        assert_eq!(trace_eval.len(), EVAL_DOM_SIZE);
    
        // commit to the trace evaluations over the larger domain using a merkle tree
        let merkle_tree = MerkleTree::<Sha3_256Backend<F>>::build(&trace_eval);
    }

}
