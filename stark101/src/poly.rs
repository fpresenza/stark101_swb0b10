use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
};
use lambdaworks_math::polynomial::{self, Polynomial};

// performs polynomial division in evaluation form.
// the obtained polynomial is the actual division if and
// only if the division remainer is zero
pub fn polynomial_division<F: IsField + IsFFTField>(
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

// performs polynomial multiplication in evaluation form.
// the obtained polynomial is the actual multiplication if
// and only if the degree of the multiplication fits in the
// domain size
pub fn polynomial_multiplication<F: IsField + IsFFTField>(
        factors: &[&Polynomial<FieldElement<F>>],
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {

    let mut product_eval = Polynomial::evaluate_offset_fft::<F>(
        factors[0], 1, Some(domain_size), offset
    ).unwrap();

    for factor in factors.iter().skip(1) {
        let evaluations = Polynomial::evaluate_offset_fft::<F>(
            factor, 1, Some(domain_size), offset
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

// performs polynomial power in evaluation form.
// the obtained polynomial is the actual power if
// and only if the degree of the power fits in the
// domain size
pub fn polynomial_power<F: IsField + IsFFTField>(
        poly: &Polynomial<FieldElement<F>>,
        power: u64,
        domain_size: usize,
        offset: &FieldElement<F>
    ) -> Polynomial<FieldElement<F>> {

    let evaluations = Polynomial::evaluate_offset_fft::<F>(
        poly, 1, Some(domain_size), offset
    ).unwrap();

    let power_eval = evaluations
            .iter()
            .map(|eval| eval.pow(power))
            .collect::<Vec<FieldElement<F>>>();

    Polynomial::interpolate_offset_fft::<F>(
        &power_eval, offset
    ).unwrap()
}

// performs polynomial folding into a new polynomial of degree
// less or equal than half the degree of the original one
pub fn fold_polynomial<F>(
    poly: &Polynomial<FieldElement<F>>,
    beta: &FieldElement<F>,
) -> Polynomial<FieldElement<F>>
where
    F: IsField,
{
    let coef = poly.coefficients();
    let even_coef: Vec<FieldElement<F>> = coef
        .iter()
        .step_by(2)
        .cloned()
        .collect();

    // odd coeficients of poly are multiplied by beta
    let odd_coef_mul_beta: Vec<FieldElement<F>> = coef
        .iter()
        .skip(1)
        .step_by(2)
        .map(|v| (v.clone()) * beta)
        .collect();

    let (even_poly, odd_poly) = polynomial::pad_with_zero_coefficients(
        &Polynomial::new(&even_coef),
        &Polynomial::new(&odd_coef_mul_beta),
    );
    even_poly + odd_poly
}