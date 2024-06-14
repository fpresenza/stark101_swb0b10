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
        // define example parameters
        let one = FE::one();
        let first_elem = one;
        let witness_elem = FE::from(3141592_u64);
        let result_elem = FE::from_hex_unchecked("6A317721EF632FF24FB815C9BBD4D4582BC7E21A43CFBDD89A8B8F0BDA68252");
        let int_dom_gen = F::get_primitive_root_of_unity(10_u64).unwrap();
        let int_dom_gen_1022 = int_dom_gen.pow(1022_u64);
        
        // create vec to hold fibonacci square sequence
        let mut fib_sq = Vec::<FE>::with_capacity(INT_DOM_SIZE);
        fib_sq.push(first_elem);
        fib_sq.push(witness_elem);
    
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
        assert_eq!(trace_poly.coefficients.len(), INT_DOM_SIZE);
        assert_eq!(trace_poly.evaluate(&one), first_elem);
        assert_eq!(trace_poly.evaluate(&int_dom_gen), witness_elem);
        assert_eq!(trace_poly.evaluate(&int_dom_gen_1022), result_elem);

        // fft-evaluate the fibonacci square sequence over a larger domain
        // of size (blow-up factor) * (interpolation domain size)
        // the offset is obtained as an outside not in the interpolation domain
        let offset = FE::from(2_u64);
        assert!(offset.pow(INT_DOM_SIZE as u64) != one);
        let trace_eval = match Polynomial::evaluate_offset_fft::<F>(
            &trace_poly, BLOWUP_FACTOR, None, &offset
        ) {
            Ok(poly) => poly,
            Err(e) => panic!("{:?}", e),
        };
        assert_eq!(trace_eval.len(), EVAL_DOM_SIZE);
    
        // commit to the trace evaluations over the larger domain using a merkle tree
        let merkle_tree = MerkleTree::<Sha3_256Backend<F>>::build(&trace_eval);

        // ===================================
        // =========|    Part 2:   |==========
        // ===== Polynomial Constraints ======
        // ===================================
        let x = Polynomial::new_monomial(one, 1);
        // initial element constraint
        let initial_constraint_poly = polynomial_division_from_evaluation(
            &trace_poly - first_elem,
            &x - one,
            Some(EVAL_DOM_SIZE),
            &offset
        );
        assert_eq!(initial_constraint_poly.coefficients.len(), 1023);

        // result element constraint
        let result_constraint_poly = polynomial_division_from_evaluation(
            &trace_poly - result_elem,
            &x - int_dom_gen_1022,
            Some(EVAL_DOM_SIZE),
            &offset
        );
        assert_eq!(result_constraint_poly.coefficients.len(), 1023);
    }

}


fn polynomial_division_from_evaluation(
        num: Polynomial<FieldElement<F>>,
        den: Polynomial<FieldElement<F>>,
        domain_size: Option<usize>,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {
    let num_eval = Polynomial::evaluate_offset_fft::<F>(
        &num, 1, domain_size, &offset
    ).unwrap();

    let den_eval = Polynomial::evaluate_offset_fft::<F>(
        &den, 1, domain_size, &offset
    ).unwrap();

    let poly_eval: Vec<_> = num_eval
        .iter()
        .zip(den_eval.iter())
        .map(|(n, d)| n / d)
        .collect();
    
    Polynomial::interpolate_fft::<F>(
        &poly_eval
    ).unwrap()

}
