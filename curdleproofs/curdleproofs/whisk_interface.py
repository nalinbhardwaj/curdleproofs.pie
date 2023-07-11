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
    PointProjective,
    point_projective_to_json,
    point_projective_from_json,
    Fr,
    BufReader,
    g1_to_bytes,
)
from py_ecc.optimized_bls12_381.optimized_curve import G1, multiply
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
from eth_typing import BLSPubkey


class WhiskTracker:
    r_G: BLSPubkey  # r * G
    k_r_G: BLSPubkey  # k * r * G

    def __init__(self, r_G: BLSPubkey, k_r_G: BLSPubkey):
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
) -> Tuple[bool, str]:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    vec_R = [pubkey_to_G1(tracker.r_G) for tracker in pre_shuffle_trackers]
    vec_S = [pubkey_to_G1(tracker.k_r_G) for tracker in pre_shuffle_trackers]

    vec_T = [pubkey_to_G1(tracker.r_G) for tracker in post_shuffle_trackers]
    vec_U = [pubkey_to_G1(tracker.k_r_G) for tracker in post_shuffle_trackers]

    ell = len(crs.vec_G)
    n_blinders = len(crs.vec_H)
    n = ell + n_blinders

    whisk_shuffle_proof = WhiskShuffleProof.from_bytes(BufReader(whisk_shuffle_proof_bytes), n)

    return whisk_shuffle_proof.proof.verify(crs, vec_R, vec_S, vec_T, vec_U, whisk_shuffle_proof.M)


def GenerateWhiskShuffleProof(
    crs: CurdleproofsCrs, pre_shuffle_trackers: Sequence[WhiskTracker]
) -> Tuple[Sequence[WhiskTracker], WhiskShuffleProofBytes]:
    permutation = list(range(len(crs.vec_G)))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [pubkey_to_G1(tracker.r_G) for tracker in pre_shuffle_trackers]
    vec_S = [pubkey_to_G1(tracker.k_r_G) for tracker in pre_shuffle_trackers]

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

    post_trackers = [WhiskTracker(G1_to_pubkey(r_G), G1_to_pubkey(k_r_G)) for r_G, k_r_G in zip(vec_T, vec_U)]

    return post_trackers, whisk_shuffle_proof.to_bytes()


SerializedWhiskTrackerProof = bytes


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSPubkey,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1_GENERATOR`.
    """
    tracker_proof_instance = TrackerOpeningProof.from_bytes(BufReader(tracker_proof))

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    return tracker_proof_instance.verify(
        transcript_verifier,
        pubkey_to_G1(tracker.k_r_G),
        pubkey_to_G1(tracker.r_G),
        pubkey_to_G1(k_commitment),
    )


def GenerateWhiskTrackerProof(
    tracker: WhiskTracker,
    k: Fr,
) -> SerializedWhiskTrackerProof:
    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=pubkey_to_G1(tracker.k_r_G),
        r_G=pubkey_to_G1(tracker.r_G),
        k_G=multiply(G1, int(k)),
        k=k,
        transcript=transcript_prover,
    )

    return opening_proof.to_bytes()
