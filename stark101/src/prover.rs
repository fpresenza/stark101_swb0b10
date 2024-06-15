use lambdaworks_math::field::{
    traits::IsFFTField,
    fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    element::FieldElement
};
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
    let g = F::get_primitive_root_of_unity(10_u64).unwrap();
    let g_to_the_1021 = g.pow(1021_u64);
    let g_to_the_1022 = g * g_to_the_1021;
    let g_to_the_1023 = g * g_to_the_1022;

    // create vec to hold fibonacci square sequence
    let mut fib_squared = Vec::<FE>::with_capacity(INT_DOM_SIZE);
    fib_squared.push(first_elem);
    fib_squared.push(witness_elem);

    for i in 2..INT_DOM_SIZE {
        let x = fib_squared[i-2];
        let y = fib_squared[i-1];
        fib_squared.push(x.pow(2_u64) + y.pow(2_u64));
    }

    // fft-interpolate the fibonacci square sequence
    let trace_poly = match Polynomial::interpolate_fft::<F>(&fib_squared) {
        Ok(poly) => poly,
        Err(e) => panic!("{:?}", e),
    };
    assert_eq!(trace_poly.coefficients.len(), INT_DOM_SIZE);
    assert_eq!(trace_poly.evaluate(&one), first_elem);
    assert_eq!(trace_poly.evaluate(&g), witness_elem);
    assert_eq!(trace_poly.evaluate(&g_to_the_1022), result_elem);

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
    let x_to_the_1024 = Polynomial::new_monomial(one, 1024);

    // initial element constraint
    let initial_constraint_poly = polynomial_division(
        &(&trace_poly - first_elem),
        &(&x - one),
        EVAL_DOM_SIZE,
        &offset
    );
    assert_eq!(initial_constraint_poly.coefficients.len(), INT_DOM_SIZE - 1);

    // result element constraint
    let result_constraint_poly = polynomial_division(
        &(&trace_poly - result_elem),
        &(&x - g_to_the_1022),
        EVAL_DOM_SIZE,
        &offset
    );
    assert_eq!(result_constraint_poly.coefficients.len(), INT_DOM_SIZE - 1);

    // trace transition constraint
    // numerator
    let trace_poly_scaled_once = trace_poly.scale(&g);
    let trace_poly_scaled_twice = trace_poly_scaled_once.scale(&g);
    let trace_poly_squared = polynomial_power(
        &trace_poly,
        2_u64,
        EVAL_DOM_SIZE,
        &offset
    );
    let trace_poly_scaled_once_squared = polynomial_power(
        &trace_poly_scaled_once,
        2_u64,
        EVAL_DOM_SIZE,
        &offset
    );
    assert!(trace_poly_squared.coefficients.len() <= 2 * INT_DOM_SIZE);
    assert!(trace_poly_scaled_once_squared.coefficients.len() <= 2 * INT_DOM_SIZE);
    assert_eq!(
        trace_poly_scaled_twice.evaluate(&g_to_the_1021),
        trace_poly_scaled_once_squared.evaluate(&g_to_the_1021) + trace_poly_squared.evaluate(&g_to_the_1021)
    );

    let numerator = polynomial_multiplication(
        &[
            &(trace_poly_scaled_twice - trace_poly_scaled_once_squared - trace_poly_squared),
            &(&x - g_to_the_1021), 
            &(&x - g_to_the_1022),
            &(&x - g_to_the_1023)
        ],
        EVAL_DOM_SIZE,
        &offset
    );

    // denominator
    let denominator = &x_to_the_1024 - one;
    let transition_constraint_poly = polynomial_division(
        &numerator,
        &denominator,
        EVAL_DOM_SIZE,
        &offset
    );
    assert!(transition_constraint_poly.coefficients.len() <= 2 * INT_DOM_SIZE);
}

// permforms polynomial division in evaluation form.
// the obtained polynomial is the actual division if and
// only if the division remainer is zero
fn polynomial_division(
        num: &Polynomial<FieldElement<F>>,
        den: &Polynomial<FieldElement<F>>,
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {

    let num_eval = Polynomial::evaluate_offset_fft::<F>(
        num, 1, Some(domain_size), offset
    ).unwrap();

    let den_eval = Polynomial::evaluate_offset_fft::<F>(
        den, 1, Some(domain_size), offset
    ).unwrap();

    let poly_eval = num_eval
        .iter()
        .zip(den_eval.iter())
        .map(|(n, d)| n / d)
        .collect::<Vec<FieldElement<F>>>();
    
    Polynomial::interpolate_offset_fft::<F>(
        &poly_eval, offset
    ).unwrap()

}

// permforms polynomial multiplication in evaluation form.
fn polynomial_multiplication(
        factors: &[&Polynomial<FieldElement<F>>],
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {

    let mut product_eval = Polynomial::evaluate_offset_fft::<F>(
        &factors[0], 1, Some(domain_size), offset
    ).unwrap();

    for i in 1..factors.len() {
        let evaluations = Polynomial::evaluate_offset_fft::<F>(
            &factors[i], 1, Some(domain_size), offset
        ).unwrap();
        product_eval = product_eval
            .iter()
            .zip(evaluations)
            .map(|(prod, eval)| prod * eval)
            .collect::<Vec<FieldElement<F>>>();
    }

    Polynomial::interpolate_offset_fft::<F>(
        &product_eval, offset
    ).unwrap()
}

// permforms polynomial power in evaluation form.
fn polynomial_power(
        poly: &Polynomial<FieldElement<F>>,
        power: u64,
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {

    let evaluations = Polynomial::evaluate_offset_fft::<F>(
        &poly, 1, Some(domain_size), offset
    ).unwrap();

    let mut power_eval = evaluations.clone();

    for _ in 1..power {
        power_eval = power_eval
            .iter()
            .zip(&evaluations)
            .map(|(pow, eval)| pow * eval)
            .collect::<Vec<FieldElement<F>>>();
    }

    Polynomial::interpolate_offset_fft::<F>(
        &power_eval, offset
    ).unwrap()
}