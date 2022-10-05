import json
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    CurdleProofsProof,
    VerifierInput,
)


def test_verifier():
    # Read fixture json
    with open("curdleproofs/fixtures/proof.json") as f:
        json_str_proof = f.read()
    with open("curdleproofs/fixtures/crs.json") as f:
        json_str_crs = f.read()
    with open("curdleproofs/fixtures/verifier_input.json") as f:
        json_str_verifier_input = f.read()

    deser_shuffle_proof = CurdleProofsProof.from_json(json_str_proof)
    deser_crs = CurdleproofsCrs.from_json(json_str_crs)
    deser_verifier_input = VerifierInput.from_json(json_str_verifier_input)

    # for i in range(50):
    #   print("iter ", i)
    verify, err = deser_shuffle_proof.verify(
        deser_crs,
        deser_verifier_input.vec_R,
        deser_verifier_input.vec_S,
        deser_verifier_input.vec_T,
        deser_verifier_input.vec_U,
        deser_verifier_input.M,
    )
    print("verify", verify)
    print("err", err)
    assert verify


def test_tracker_opening_proof():
    G = G1
    k = generate_blinders(1)[0]
    r = generate_blinders(1)[0]

    k_G = multiply(G, int(k))
    r_G = multiply(G, int(r))
    k_r_G = multiply(r_G, int(k))

    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=k_r_G, r_G=r_G, k_G=k_G, G=G, k=k, transcript=transcript_prover
    )

    json_str_proof = opening_proof.to_json()
    print("json_str_proof", json_str_proof)

    deser_proof = TrackerOpeningProof.from_json(json_str_proof)

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    assert deser_proof.verify(transcript_verifier)
