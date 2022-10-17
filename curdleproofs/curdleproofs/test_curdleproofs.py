import json
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    CurdleProofsProof,
    VerifierInput,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from py_ecc.optimized_bls12_381.optimized_curve import G1


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
    # Read fixture json
    with open("curdleproofs/fixtures/tracker_opening_proof.json") as f:
        json_str_proof = f.read()
    deser_proof = TrackerOpeningProof.from_json(json_str_proof)

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    assert deser_proof.verify(transcript_verifier)
