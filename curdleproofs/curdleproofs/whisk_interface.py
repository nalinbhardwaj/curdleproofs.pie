import random
import json
from typing import Container, Sequence, Tuple, Type, TypeVar
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    N_BLINDERS,
    CurdleProofsProof,
    shuffle_permute_and_commit_input,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import (
    PointAffine as BLSG1Point,
    PointProjective,
    affine_to_projective,
    point_projective_to_json,
    point_projective_from_json,
    Fr,
)
from py_ecc.optimized_bls12_381.optimized_curve import G1, normalize, multiply


class WhiskTracker:
    r_G: PointProjective  # r * G
    k_r_G: PointProjective  # k * r * G

    def __init__(self, r_G: PointProjective, k_r_G: PointProjective):
        self.r_G = r_G
        self.k_r_G = k_r_G

T_WhiskShuffleProof = TypeVar("T_WhiskShuffleProof", bound="WhiskShuffleProof")

class WhiskShuffleProof:
    M: PointProjective
    proof: CurdleProofsProof

    def __init__(self, M: PointProjective, proof: CurdleProofsProof):
        self.M = M
        self.proof = proof

    def to_json(self):
        return {
            "M": point_projective_to_json(self.M),
            "proof": self.proof.to_json(),
        }

    @classmethod
    def from_json(cls: Type[T_WhiskShuffleProof], data) -> T_WhiskShuffleProof:
        return cls(
            M=point_projective_from_json(data["M"]),
            proof=CurdleProofsProof.from_json(data["proof"]),
        )


WhiskShuffleProofBytes = bytes


def IsValidWhiskShuffleProof(
    crs: CurdleproofsCrs,
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    whisk_shuffle_proof_bytes: WhiskShuffleProofBytes,
) -> Tuple[bool, str]:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    vec_R = [tracker.r_G for tracker in pre_shuffle_trackers]
    vec_S = [tracker.k_r_G for tracker in pre_shuffle_trackers]

    vec_T = [tracker.r_G for tracker in post_shuffle_trackers]
    vec_U = [tracker.k_r_G for tracker in post_shuffle_trackers]

    whisk_shuffle_proof = WhiskShuffleProof.from_json(json.loads(whisk_shuffle_proof_bytes.decode()))

    return whisk_shuffle_proof.proof.verify(crs, vec_R, vec_S, vec_T, vec_U, whisk_shuffle_proof.M)


def GenerateWhiskShuffleProof(
    crs: CurdleproofsCrs, pre_shuffle_trackers: Sequence[WhiskTracker]
) -> Tuple[Sequence[WhiskTracker], WhiskShuffleProofBytes]:
    permutation = list(range(len(crs.vec_G)))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [tracker.r_G for tracker in pre_shuffle_trackers]
    vec_S = [tracker.k_r_G for tracker in pre_shuffle_trackers]

    vec_T, vec_U, M, vec_m_blinders = shuffle_permute_and_commit_input(
        crs, vec_R, vec_S, permutation, k
    )

    shuffle_proof: CurdleProofsProof = CurdleProofsProof.new(
        crs=crs,
        vec_R=vec_R,
        vec_S=vec_S,
        vec_T=vec_T,
        vec_U=vec_U,
        M=M,
        permutation=permutation,
        k=k,
        vec_m_blinders=vec_m_blinders,
    )
    whisk_shuffle_proof = WhiskShuffleProof(M, shuffle_proof)

    post_trackers = [WhiskTracker(r_G, k_r_G) for r_G, k_r_G in zip(vec_T, vec_U)]

    return post_trackers, json.dumps(whisk_shuffle_proof.to_json()).encode()


SerializedWhiskTrackerProof = bytes


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: PointProjective,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1_GENERATOR`.
    """
    tracker_proof_instance = TrackerOpeningProof.from_json(json.loads(tracker_proof.decode()))

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    return tracker_proof_instance.verify(
        transcript_verifier,
        tracker.k_r_G,
        tracker.r_G,
        k_commitment,
    )


def GenerateWhiskTrackerProof(
    tracker: WhiskTracker,
    k: Fr,
) -> SerializedWhiskTrackerProof:
    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=tracker.k_r_G,
        r_G=tracker.r_G,
        k_G=multiply(G1, int(k)),
        G=G1,
        k=k,
        transcript=transcript_prover,
    )

    return json.dumps(opening_proof.to_json()).encode()
