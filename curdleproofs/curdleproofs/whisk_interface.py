import random
from typing import Sequence, Tuple, Type, TypeVar
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    CurdleProofsProof,
    shuffle_permute_and_commit_input,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import (
    point_projective_to_json,
    point_projective_from_json,
    point_projective_to_bytes,
    point_projective_from_bytes,
    BufReader,
    g1_to_bytes,
    random_scalar,
    G1,
    BLSPubkey
)
from py_arkworks_bls12381 import G1Point, Scalar


class WhiskTracker:
    r_G: BLSPubkey  # r * G
    k_r_G: BLSPubkey  # k * r * G

    def __init__(self, r_G: BLSPubkey, k_r_G: BLSPubkey):
        self.r_G = r_G
        self.k_r_G = k_r_G


T_WhiskShuffleProof = TypeVar("T_WhiskShuffleProof", bound="WhiskShuffleProof")


class WhiskShuffleProof:
    M: G1Point
    proof: CurdleProofsProof

    def __init__(self, M: G1Point, proof: CurdleProofsProof):
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

    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.M),
            self.proof.to_bytes(),
        ])

    @classmethod
    def from_bytes(cls: Type[T_WhiskShuffleProof], b: BufReader, n: int) -> T_WhiskShuffleProof:
        return cls(
            M=b.read_g1(),
            proof=CurdleProofsProof.from_bytes(b, n),
        )


WhiskShuffleProofBytes = bytes


def IsValidWhiskShuffleProof(
    crs: CurdleproofsCrs,
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    whisk_shuffle_proof_bytes: WhiskShuffleProofBytes,
) -> bool:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    try:
        AssertIsValidWhiskShuffleProof(crs, pre_shuffle_trackers, post_shuffle_trackers, whisk_shuffle_proof_bytes)
        return True
    except:  # noqa: E722
        return False


def AssertIsValidWhiskShuffleProof(
    crs: CurdleproofsCrs,
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    whisk_shuffle_proof_bytes: WhiskShuffleProofBytes,
):
    vec_R = [point_projective_from_bytes(tracker.r_G) for tracker in pre_shuffle_trackers]
    vec_S = [point_projective_from_bytes(tracker.k_r_G) for tracker in pre_shuffle_trackers]

    vec_T = [point_projective_from_bytes(tracker.r_G) for tracker in post_shuffle_trackers]
    vec_U = [point_projective_from_bytes(tracker.k_r_G) for tracker in post_shuffle_trackers]

    ell = len(crs.vec_G)
    n_blinders = len(crs.vec_H)
    n = ell + n_blinders

    whisk_shuffle_proof = WhiskShuffleProof.from_bytes(BufReader(whisk_shuffle_proof_bytes), n)

    whisk_shuffle_proof.proof.verify(crs, vec_R, vec_S, vec_T, vec_U, whisk_shuffle_proof.M)


def GenerateWhiskShuffleProof(
    crs: CurdleproofsCrs, pre_shuffle_trackers: Sequence[WhiskTracker]
) -> Tuple[Sequence[WhiskTracker], WhiskShuffleProofBytes]:
    permutation = list(range(len(crs.vec_G)))
    random.shuffle(permutation)
    k = random_scalar()

    vec_R = [point_projective_from_bytes(tracker.r_G) for tracker in pre_shuffle_trackers]
    vec_S = [point_projective_from_bytes(tracker.k_r_G) for tracker in pre_shuffle_trackers]

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

    post_trackers = [WhiskTracker(BLSPubkey(point_projective_to_bytes(r_G)), BLSPubkey(point_projective_to_bytes(k_r_G))) for r_G, k_r_G in zip(vec_T, vec_U)]

    return post_trackers, whisk_shuffle_proof.to_bytes()


SerializedWhiskTrackerProof = bytes


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSPubkey,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1`.
    """
    try:
        AssertIsValidWhiskOpeningProof(tracker, k_commitment, tracker_proof)
        return True
    except:  # noqa: E722
        return False


def AssertIsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSPubkey,
    tracker_proof: SerializedWhiskTrackerProof,
):
    tracker_proof_instance = TrackerOpeningProof.from_bytes(BufReader(tracker_proof))

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    tracker_proof_instance.verify(
        transcript_verifier,
        point_projective_from_bytes(tracker.k_r_G),
        point_projective_from_bytes(tracker.r_G),
        point_projective_from_bytes(k_commitment),
    )


def GenerateWhiskTrackerProof(
    tracker: WhiskTracker,
    k: Scalar,
) -> SerializedWhiskTrackerProof:
    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=point_projective_from_bytes(tracker.k_r_G),
        r_G=point_projective_from_bytes(tracker.r_G),
        k_G=G1 * k,
        k=k,
        transcript=transcript_prover,
    )

    return opening_proof.to_bytes()
