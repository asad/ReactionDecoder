package org.openscience.smsd.algorithm.vflib.vf2.mcs;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Given a (subgraph-)isomorphism state this class can lazily iterate over the
 * mappings in a non-recursive manner. The class currently implements and {@link
 * Iterator} but is better suited to the {@code Stream} class (which will be
 * available in JDK 8).
 *
 * @author John May
 * @author Syed Asad Rahman
 * 
 */
final class StateStream implements Iterator<Collection<int[]>> {

    /**
     * A mapping state.
     */
    private final State state;

    /**
     * The stack replaces the call-stack in a recursive matcher.
     */
    private final CandidateStack stack;

    /**
     * Current candidates.
     */
    private int n = 0, m = -1;

    /**
     * The next mapping.
     */
    private Collection<int[]> next;

    /**
     * Create a stream for the provided state.
     *
     * @param state the state to stream over
     */
    StateStream(final State state) {
        this.state = state;
        this.stack = new CandidateStack(state.maxQueryCandidate());
        this.next = state.maxQueryCandidate() == 0 || state.maxTargetCandidate() == 0 ? null : findNext(); // first-mapping
    }

    /**
     * @inheritDoc
     */
    @Override
    public boolean hasNext() {
        return next != null;
    }

    /**
     * @inheritDoc
     */
    @Override
    public Collection<int[]> next() {
        Collection<int[]> ret = next;
        next = findNext();
        return ret;
    }

    /**
     * @inheritDoc
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("a graph matching cannot be removed");
    }

    /**
     * Finds the next mapping from the current state.
     *
     * @return the next state (or null if none)
     */
    private Collection<int[]> findNext() {

        List<int[]> types = new LinkedList<>();
        /*
             * Return maximum match
         */
        while (map()) {
            if (state.size() <= state.maxQueryCandidate() && state.size() > 0) {
                if (types.isEmpty()) {
                    types.add(state.mapping());
                } else if (matchLength(types.iterator().next()) < matchLength(state.mapping())) {
                    types.clear();
                    types.add(state.mapping());
                } else if (matchLength(types.iterator().next()) == matchLength(state.mapping())) {
                    if (!solutionPresent(types, state.mapping())) {
                        types.add(state.mapping());
                    }
                }
            }
        }
//        System.out.println("\n\nSolutions from VF");
//        types.stream().forEach((s) -> {
//            System.out.println("Solutions " + Arrays.toString(s));
//        });
        return types.isEmpty() ? null : types;

    }

    private boolean solutionPresent(Collection<int[]> storedSolutions, int[] potentialSolution) {
        for (int[] solution : storedSolutions) {
            boolean equalsFlag = Arrays.equals(solution, potentialSolution);
            if (equalsFlag) {
                return equalsFlag;
            }
        }
        return false;
    }

    private int matchLength(int[] match) {
        int counter = 0;
        for (int i : match) {
            if (i >= 0) {
                counter++;
            }
        }
        return counter;
    }

    /**
     * Progress the state-machine - the function return false when a mapping is
     * found on the mapping is done.
     *
     * @return the state is partial
     */
    private boolean map() {

        // backtrack - we've tried all possible n or m, backTrack the last mapping
        if ((n == state.maxQueryCandidate() || m == state.maxTargetCandidate()) && !stack.empty()) {
            state.backTrack(n = stack.popN(), m = stack.popM());
        }

        while ((m = state.nextCandidate(n, m)) < state.maxTargetCandidate()) {
            if (state.addMapping(n, m)) {
                stack.push(n, m);
                n = state.hasNextCandidate(-1);
                m = -1;
                return n < state.maxQueryCandidate();
            }
        }

        return state.size() > 0 || m < state.maxTargetCandidate();
    }

    /**
     * A fixed size stack to keep track of which vertices are mapped. This stack
     * allows us to turn the recursive algorithms it to lazy iterating mappers.
     * A reclusive call is usually implemented as call-stack which stores the
     * variable in each subroutine invocation. For the mapping we actually only
     * need store the candidates.
     */
    private final class CandidateStack {

        /**
         * Candidate storage.
         */
        private final int[] ns, ms;

        /**
         * Size of each stack.
         */
        private int nSize, mSize;

        private CandidateStack(int capacity) {
            ns = new int[capacity];
            ms = new int[capacity];
        }

        /**
         * Push a candidate mapping on to the stack.
         *
         * @param n query candidate
         * @param m target candidate
         */
        void push(int n, int m) {
            ns[nSize++] = n;
            ms[mSize++] = m;
        }

        /**
         * Pops the G1 candidate.
         *
         * @return the previous 'n' candidate
         */
        int popN() {
            return ns[--nSize];
        }

        /**
         * Pops the G2 candidate.
         *
         * @return the previous 'm' candidate
         */
        int popM() {
            return ms[--mSize];
        }

        /**
         * Is the stack empty - if so no candidates can be popped.
         *
         * @return
         */
        boolean empty() {
            return nSize == 0 && mSize == 0;
        }
    }
}
