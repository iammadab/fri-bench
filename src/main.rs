use std::time::Instant;
use plonky2::fri::oracle::PolynomialBatch;
use plonky2::fri::reduction_strategies::FriReductionStrategy;
use plonky2::fri::structure::{
    FriBatchInfo, FriInstanceInfo, FriOpeningBatch, FriOpenings, FriOracleInfo, FriPolynomialInfo,
};
use plonky2::fri::verifier::verify_fri_proof;
use plonky2::fri::FriConfig;
use plonky2::iop::challenger::Challenger;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
use plonky2::util::timing::TimingTree;
use plonky2_field::goldilocks_field::GoldilocksField;
use plonky2_field::polynomial::PolynomialValues;
use plonky2_field::types::Sample;

fn bench_commit_and_prove(n_vars: usize) {
    let mut timing = TimingTree::default();
    let degree_bits = n_vars;
    let hiding = false;
    let evaluations = GoldilocksField::rand_vec(1 << degree_bits);

    let fri_config = FriConfig {
        rate_bits: 1,
        cap_height: 4,
        proof_of_work_bits: 16,
        reduction_strategy: FriReductionStrategy::ConstantArityBits(4, 5),
        num_query_rounds: 84,
    };
    // let fri_config = FriConfig {
    //         rate_bits: 1,
    //         cap_height: 5,
    //         proof_of_work_bits: 0,
    //         reduction_strategy: FriReductionStrategy::Fixed(vec![1, 2, 1]),
    //         num_query_rounds: 10,
    // };
    let fri_params = fri_config.fri_params(degree_bits, hiding);

    let polynomial_values = PolynomialValues::new(evaluations);

    let now = Instant::now();
    // TODO: get the rate bit from the fri config
    // TODO: maybe benchmark also for blinding (structure should be generic enough to easily allow for this)
    // TODO: what is cap height
    // TODO: any problem with setting fft tables to none?
    // TODO: figure out the difference when you tweak the cap height
    let polynomial_batch: PolynomialBatch<GoldilocksField, PoseidonGoldilocksConfig, 2> =
        PolynomialBatch::from_values(
            vec![polynomial_values],
            fri_config.rate_bits,
            hiding,
            fri_config.cap_height,
            &mut timing,
            None,
        );
    println!("commit-time: {:?}ms", now.elapsed().as_millis());

    let mut challenger =
        Challenger::<GoldilocksField, <PoseidonGoldilocksConfig as GenericConfig<2>>::Hasher>::new(
        );
    challenger.observe_cap(&polynomial_batch.merkle_tree.cap);

    // challenger.get_n_challenges(2);
    let challenge_point = challenger.get_extension_challenge::<2>();

    // the polynomial is opened at some random point
    let point_batch: FriBatchInfo<GoldilocksField, 2> = FriBatchInfo {
        point: challenge_point,
        polynomials: vec![FriPolynomialInfo {
            oracle_index: 0,
            polynomial_index: 0,
        }],
    };

    let fri_instance = FriInstanceInfo {
        oracles: vec![FriOracleInfo {
            num_polys: 1,
            blinding: hiding,
        }],
        batches: vec![point_batch],
    };

    let now = Instant::now();
    let proof = PolynomialBatch::<GoldilocksField, PoseidonGoldilocksConfig, 2>::prove_openings(
        &fri_instance,
        &[&polynomial_batch],
        &mut challenger,
        &fri_params,
        &mut timing,
    );
    println!("prove-time: {:?}ms", now.elapsed().as_millis());

    // how do I verify
    // seems I need to construct the fri opening batch myself
    // one polynomial evaluated at the challenge point
    let fri_opening: FriOpenings<GoldilocksField, 2> = FriOpenings {
        batches: vec![FriOpeningBatch {
            values: vec![polynomial_batch.polynomials[0]
                .to_extension::<2>()
                .eval(challenge_point)],
        }],
    };

    let fri_challenges = challenger.fri_challenges::<PoseidonGoldilocksConfig, 2>(
        &proof.commit_phase_merkle_caps,
        &proof.final_poly,
        proof.pow_witness,
        degree_bits,
        &fri_config,
    );

    let now = Instant::now();
    // merkle caps seem to be a collection of caps
    let verification_result = verify_fri_proof::<GoldilocksField, PoseidonGoldilocksConfig, 2>(
        &fri_instance,
        &fri_opening,
        &fri_challenges,
        &[polynomial_batch.merkle_tree.cap],
        &proof,
        &fri_params,
    );
    println!("verification-time: {:?}ms", now.elapsed().as_millis());

    dbg!(verification_result);
}

fn bench_commit_and_prove_60_20_var_poly() {
    let mut timing = TimingTree::default();
    let degree_bits = 20;
    let hiding = false;

    let fri_config = FriConfig {
        rate_bits: 1,
        cap_height: 4,
        proof_of_work_bits: 16,
        reduction_strategy: FriReductionStrategy::ConstantArityBits(4, 5),
        num_query_rounds: 84,
    };
    let fri_params = fri_config.fri_params(degree_bits, hiding);

    let polynomial_values = (0..60).map(|_| PolynomialValues::new(GoldilocksField::rand_vec(1 << degree_bits))).collect::<Vec<_>>();

    let now = Instant::now();
    let polynomial_batch: PolynomialBatch<GoldilocksField, PoseidonGoldilocksConfig, 2> = PolynomialBatch::from_values(
        polynomial_values,
        fri_config.rate_bits,
        hiding,
        fri_config.cap_height,
        &mut timing,
        None
    );
    println!("commit-time: {:?}ms", now.elapsed().as_millis());

    let mut challenger =
        Challenger::<GoldilocksField, <PoseidonGoldilocksConfig as GenericConfig<2>>::Hasher>::new(
        );

    let challenge_point = challenger.get_extension_challenge::<2>();

    let point_batch: FriBatchInfo<GoldilocksField, 2> = FriBatchInfo {
        point: challenge_point,
        polynomials: FriPolynomialInfo::from_range(0, 0..60)
    };

    let fri_instance = FriInstanceInfo {
        oracles: vec![FriOracleInfo {
            num_polys: 60,
            blinding: hiding
        }],
        batches: vec![point_batch]
    };

    let now = Instant::now();
    let proof = PolynomialBatch::<GoldilocksField, PoseidonGoldilocksConfig, 2>::prove_openings(
        &fri_instance,
        &[&polynomial_batch],
        &mut challenger,
        &fri_params,
        &mut timing
    );
    println!("prove-time: {:?}ms", now.elapsed().as_millis());
}

fn main() {
    // commit - prove - verify cycle for n_vars = 17..=22
    for n in 17..=22 {
        println!("benchmarks for n_vars = {n}");
        bench_commit_and_prove(n);
        println!("");
    }

    println!("benchmarks for 60 2^20 evaluations");
    bench_commit_and_prove_60_20_var_poly();
}
