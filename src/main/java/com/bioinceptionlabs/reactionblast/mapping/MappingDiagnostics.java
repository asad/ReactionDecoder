/*
 * Copyright (C) 2003-2026 Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */
package com.bioinceptionlabs.reactionblast.mapping;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Per-reaction diagnostics for mapping hot-path benchmarking.
 */
public final class MappingDiagnostics {

    private static final ConcurrentMap<String, ReactionStats> REACTIONS = new ConcurrentHashMap<>();

    private MappingDiagnostics() {
    }

    public static void resetReaction(String reactionId) {
        if (reactionId != null) {
            REACTIONS.remove(reactionId);
        }
    }

    public static int recordMatcherInvocation(String reactionId, String algorithm,
            long candidatePairs, long uniquePairs,
            long identitySkips, long ratioSkips, long tanimotoSkips,
            long scheduledJobs) {
        return reactionStats(reactionId)
                .algorithmStats(algorithm)
                .recordMatcherInvocation(candidatePairs, uniquePairs,
                        identitySkips, ratioSkips, tanimotoSkips, scheduledJobs);
    }

    public static void recordMatcherCompletion(String reactionId, String algorithm,
            int invocationIndex, long replayedMappings, long elapsedMillis) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .invocation(invocationIndex)
                .complete(replayedMappings, elapsedMillis);
    }

    public static void recordSubstructureSearch(String reactionId, String algorithm, int invocationIndex) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .invocation(invocationIndex)
                .substructureSearches.incrementAndGet();
    }

    public static void recordMcsCacheHit(String reactionId, String algorithm, int invocationIndex) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .invocation(invocationIndex)
                .cacheHits.incrementAndGet();
    }

    public static void recordActualMcsSearch(String reactionId, String algorithm, int invocationIndex) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .invocation(invocationIndex)
                .actualMcsSearches.incrementAndGet();
    }

    public static void recordQuickMappingCall(String reactionId, String algorithm) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .quickMappingCalls.incrementAndGet();
    }

    public static void recordQuickMappingCacheHit(String reactionId, String algorithm) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .quickMappingCacheHits.incrementAndGet();
    }

    public static void recordQuickMappingSearch(String reactionId, String algorithm) {
        reactionStats(reactionId)
                .algorithmStats(algorithm)
                .quickMappingSearches.incrementAndGet();
    }

    public static void recordMappingPhase(String reactionId, long elapsedMillis) {
        reactionStats(reactionId).mappingPhaseMillis.set(elapsedMillis);
    }

    public static void recordEvaluationPhase(String reactionId, long elapsedMillis) {
        reactionStats(reactionId).evaluationPhaseMillis.set(elapsedMillis);
    }

    public static ReactionSnapshot snapshot(String reactionId) {
        ReactionStats stats = REACTIONS.get(reactionId);
        return stats == null
                ? new ReactionSnapshot(reactionId, 0L, 0L, Collections.emptyList())
                : stats.snapshot(reactionId);
    }

    private static ReactionStats reactionStats(String reactionId) {
        return REACTIONS.computeIfAbsent(
                reactionId == null ? "UNKNOWN_REACTION" : reactionId,
                ignored -> new ReactionStats());
    }

    private static final class ReactionStats {

        private final ConcurrentMap<String, AlgorithmStats> algorithms = new ConcurrentHashMap<>();
        private final AtomicLong mappingPhaseMillis = new AtomicLong();
        private final AtomicLong evaluationPhaseMillis = new AtomicLong();

        private AlgorithmStats algorithmStats(String algorithm) {
            String key = algorithm == null ? "UNKNOWN" : algorithm;
            return algorithms.computeIfAbsent(key, AlgorithmStats::new);
        }

        private ReactionSnapshot snapshot(String reactionId) {
            List<AlgorithmSnapshot> algorithmSnapshots = new ArrayList<>();
            for (AlgorithmStats stats : algorithms.values()) {
                algorithmSnapshots.add(stats.snapshot());
            }
            algorithmSnapshots.sort(Comparator.comparing(snapshot -> snapshot.algorithm));
            return new ReactionSnapshot(
                    reactionId,
                    mappingPhaseMillis.get(),
                    evaluationPhaseMillis.get(),
                    algorithmSnapshots);
        }
    }

    private static final class AlgorithmStats {

        private final String algorithm;
        private final AtomicInteger matcherInvocationCounter = new AtomicInteger();
        private final ConcurrentMap<Integer, InvocationStats> invocations = new ConcurrentHashMap<>();
        private final AtomicLong quickMappingCalls = new AtomicLong();
        private final AtomicLong quickMappingCacheHits = new AtomicLong();
        private final AtomicLong quickMappingSearches = new AtomicLong();

        private AlgorithmStats(String algorithm) {
            this.algorithm = algorithm;
        }

        private int recordMatcherInvocation(long candidatePairs, long uniquePairs,
                long identitySkips, long ratioSkips, long tanimotoSkips,
                long scheduledJobs) {
            int invocationIndex = matcherInvocationCounter.incrementAndGet();
            InvocationStats invocation = invocation(invocationIndex);
            invocation.candidatePairs.set(candidatePairs);
            invocation.uniquePairs.set(uniquePairs);
            invocation.identitySkips.set(identitySkips);
            invocation.ratioSkips.set(ratioSkips);
            invocation.tanimotoSkips.set(tanimotoSkips);
            invocation.scheduledJobs.set(scheduledJobs);
            return invocationIndex;
        }

        private InvocationStats invocation(int invocationIndex) {
            return invocations.computeIfAbsent(invocationIndex, InvocationStats::new);
        }

        private AlgorithmSnapshot snapshot() {
            List<MatcherInvocationSnapshot> invocationSnapshots = new ArrayList<>();
            for (InvocationStats stats : invocations.values()) {
                invocationSnapshots.add(stats.snapshot());
            }
            invocationSnapshots.sort(Comparator.comparingInt(snapshot -> snapshot.invocationIndex));
            return new AlgorithmSnapshot(
                    algorithm,
                    quickMappingCalls.get(),
                    quickMappingCacheHits.get(),
                    quickMappingSearches.get(),
                    invocationSnapshots);
        }
    }

    private static final class InvocationStats {

        private final int invocationIndex;
        private final AtomicLong candidatePairs = new AtomicLong();
        private final AtomicLong uniquePairs = new AtomicLong();
        private final AtomicLong identitySkips = new AtomicLong();
        private final AtomicLong ratioSkips = new AtomicLong();
        private final AtomicLong tanimotoSkips = new AtomicLong();
        private final AtomicLong scheduledJobs = new AtomicLong();
        private final AtomicLong replayedMappings = new AtomicLong();
        private final AtomicLong substructureSearches = new AtomicLong();
        private final AtomicLong cacheHits = new AtomicLong();
        private final AtomicLong actualMcsSearches = new AtomicLong();
        private final AtomicLong elapsedMillis = new AtomicLong();

        private InvocationStats(int invocationIndex) {
            this.invocationIndex = invocationIndex;
        }

        private void complete(long replayedMappings, long elapsedMillis) {
            this.replayedMappings.set(replayedMappings);
            this.elapsedMillis.set(elapsedMillis);
        }

        private MatcherInvocationSnapshot snapshot() {
            return new MatcherInvocationSnapshot(
                    invocationIndex,
                    candidatePairs.get(),
                    uniquePairs.get(),
                    identitySkips.get(),
                    ratioSkips.get(),
                    tanimotoSkips.get(),
                    scheduledJobs.get(),
                    replayedMappings.get(),
                    substructureSearches.get(),
                    cacheHits.get(),
                    actualMcsSearches.get(),
                    elapsedMillis.get());
        }
    }

    public static final class ReactionSnapshot {

        public final String reactionId;
        public final long mappingPhaseMillis;
        public final long evaluationPhaseMillis;
        public final List<AlgorithmSnapshot> algorithms;

        public ReactionSnapshot(String reactionId,
                long mappingPhaseMillis,
                long evaluationPhaseMillis,
                List<AlgorithmSnapshot> algorithms) {
            this.reactionId = reactionId;
            this.mappingPhaseMillis = mappingPhaseMillis;
            this.evaluationPhaseMillis = evaluationPhaseMillis;
            this.algorithms = Collections.unmodifiableList(new ArrayList<>(algorithms));
        }
    }

    public static final class AlgorithmSnapshot {

        public final String algorithm;
        public final long quickMappingCalls;
        public final long quickMappingCacheHits;
        public final long quickMappingSearches;
        public final List<MatcherInvocationSnapshot> invocations;

        public AlgorithmSnapshot(String algorithm,
                long quickMappingCalls,
                long quickMappingCacheHits,
                long quickMappingSearches,
                List<MatcherInvocationSnapshot> invocations) {
            this.algorithm = algorithm;
            this.quickMappingCalls = quickMappingCalls;
            this.quickMappingCacheHits = quickMappingCacheHits;
            this.quickMappingSearches = quickMappingSearches;
            this.invocations = Collections.unmodifiableList(new ArrayList<>(invocations));
        }
    }

    public static final class MatcherInvocationSnapshot {

        public final int invocationIndex;
        public final long candidatePairs;
        public final long uniquePairs;
        public final long identitySkips;
        public final long ratioSkips;
        public final long tanimotoSkips;
        public final long scheduledJobs;
        public final long replayedMappings;
        public final long substructureSearches;
        public final long cacheHits;
        public final long actualMcsSearches;
        public final long elapsedMillis;

        public MatcherInvocationSnapshot(int invocationIndex,
                long candidatePairs,
                long uniquePairs,
                long identitySkips,
                long ratioSkips,
                long tanimotoSkips,
                long scheduledJobs,
                long replayedMappings,
                long substructureSearches,
                long cacheHits,
                long actualMcsSearches,
                long elapsedMillis) {
            this.invocationIndex = invocationIndex;
            this.candidatePairs = candidatePairs;
            this.uniquePairs = uniquePairs;
            this.identitySkips = identitySkips;
            this.ratioSkips = ratioSkips;
            this.tanimotoSkips = tanimotoSkips;
            this.scheduledJobs = scheduledJobs;
            this.replayedMappings = replayedMappings;
            this.substructureSearches = substructureSearches;
            this.cacheHits = cacheHits;
            this.actualMcsSearches = actualMcsSearches;
            this.elapsedMillis = elapsedMillis;
        }
    }
}
