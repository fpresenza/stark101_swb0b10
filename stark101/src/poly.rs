use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
};
use lambdaworks_math::polynomial::Polynomial;

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
        &poly, 1, Some(domain_size), offset
    ).unwrap();

    let power_eval = evaluations
            .iter()
            .map(|eval| eval.pow(power))
            .collect::<Vec<FieldElement<F>>>();

    Polynomial::interpolate_offset_fft::<F>(
        &power_eval, offset
    ).unwrap()
}