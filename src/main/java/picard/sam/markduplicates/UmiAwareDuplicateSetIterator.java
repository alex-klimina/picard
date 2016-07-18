/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import picard.PicardException;

import java.util.*;

/**
 * Created by fleharty on 5/23/16.
 */
public class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private int editDistanceToJoin;
    private boolean addInferredUmi;
    private String umiTag;
    private String inferredUmiTag;

    public UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int editDistanceToJoin, final boolean addInferredUmi, final String umiTag, final String inferredUmiTag) {
        this.wrappedIterator = wrappedIterator;
        this.editDistanceToJoin = editDistanceToJoin;
        this.addInferredUmi = addInferredUmi;
        this.umiTag = umiTag;
        this.inferredUmiTag = inferredUmiTag;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        wrappedIterator.close();
    }

    @Override
    public boolean hasNext() {
        return nextSetsIterator.hasNext() || wrappedIterator.hasNext();
    }

    @Override
    public DuplicateSet next() {
        if(!nextSetsIterator.hasNext())
            process(wrappedIterator.next());

        return nextSetsIterator.next();
    }

    // Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
    // and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
    private void process(final DuplicateSet set) {

        List<SAMRecord> records = set.getRecords();

        // If any records are missing the UMI_TAG proceed as if there were no UMIs
        // and return nextSetsIterator without breaking it up into smaller sets.
        for(SAMRecord rec : records) {
            if(rec.getStringAttribute(umiTag) == null) {
                nextSetsIterator = Collections.singleton(set).iterator();
                return;
            }
        }

        // Sort records by RX tag
        Collections.sort(records, new Comparator<SAMRecord>() {
            @Override
            public int compare(final SAMRecord lhs, final SAMRecord rhs) {
                return (lhs.getStringAttribute(umiTag)).compareTo(rhs.getStringAttribute(umiTag));
            }
        });

        // Locate records that have identical UMI sequences
        Set<String> uniqueObservedUMIs = new HashSet();
        for(SAMRecord rec : records) {
            uniqueObservedUMIs.add(rec.getStringAttribute(umiTag));
        }

        List<String> observedUmis = new ArrayList<>(uniqueObservedUMIs);
        List<List<Integer>> adjacencyList = new ArrayList<>();
        List<Integer> duplicateSets = new ArrayList<>();

        // Construct adjacencyList and groups

        for(String umi1 : observedUmis) {
            int adjacencyListIndex = observedUmis.indexOf(umi1);
            adjacencyList.add(adjacencyListIndex, new ArrayList<>());
            duplicateSets.add(adjacencyListIndex, 0);

            for(String umi2 : observedUmis) {
                int adjacencyListIndex2 = observedUmis.indexOf(umi2);
                if( getEditDistance(observedUmis.get(adjacencyListIndex), observedUmis.get(adjacencyListIndex2)) <= editDistanceToJoin) {
                    adjacencyList.get(adjacencyListIndex).add(adjacencyListIndex2);
                }
            }
        }

        // Join Groups
        int nDuplicateSets = 0;
        for(int i = 0;i < adjacencyList.size();i++) {
            // Have I visited this yet?
            if(duplicateSets.get(i) == 0) {
                // No, I haven't yet seen this
                nDuplicateSets++; // We've now seen a new group

                // Depth first search on adjacencyList, setting all the values to a group
                mergeDuplicateSets(i, nDuplicateSets, duplicateSets, adjacencyList);
            }
        }

        // Construct DuplicateSetList
        List<DuplicateSet> duplicateSetList = new ArrayList<>();
        for(int i = 0;i < nDuplicateSets;i++) {
            DuplicateSet e = new DuplicateSet();
            duplicateSetList.add(e);
        }

        // Assign each record to a duplicate set
        for(SAMRecord rec : records) {
            String umi = records.get(records.indexOf(rec)).getStringAttribute(umiTag);

            // Figure out which group this belongs to
            int recordGroup = duplicateSets.get(observedUmis.indexOf(umi));
            duplicateSetList.get(recordGroup-1).add(records.get(records.indexOf(rec)));

        }

        // Optionally add the inferred Umi
        if(addInferredUmi) {
            for(DuplicateSet ds : duplicateSetList) {
                // For each duplicate set identify the most common UMI
                List<SAMRecord> recordsFromDuplicateSet = ds.getRecords();

                // Count the number of times each UMI appears
                Map<String, Long> umiCounts = new HashMap<>();
                for(SAMRecord record : recordsFromDuplicateSet) {
                    String currentUmi = record.getStringAttribute(umiTag);
                    Long currentCount = umiCounts.get(currentUmi);
                    if(currentCount == null) {
                        umiCounts.put(currentUmi, new Long(1));
                    } else {
                        umiCounts.put(currentUmi, currentCount + 1);
                    }
                }

                // Find most common UMI
                Map.Entry<String, Long> maxEntry = null;
                for(Map.Entry<String, Long> entry : umiCounts.entrySet()) {
                    if(maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0) {
                        maxEntry = entry;
                    }
                }

                // Assign the inferred UMI to all reads in the current group
                for(SAMRecord record : recordsFromDuplicateSet) {
                    record.setAttribute(inferredUmiTag, maxEntry.getKey());
                }
            }
        }

        nextSetsIterator = duplicateSetList.iterator();
    }

    // Assigns all connected components in the adjacencyList to the same duplicate set.
    // Each duplicate set will be labeled as a unique integer.
    private void mergeDuplicateSets(final int index, final int duplicateSet, final List<Integer> duplicateSets, final List<List<Integer>> adjacencyList) {
        // A new group is created when the group index is 0 indicating it hasn't been seen before.
        // If it has been seen before, it won't be 0, and we will skip it.
        // This is guaranteed to terminate since we traverse each member of the adjacency list once and only once,
        // if we see it again we do not recurse.
        if(duplicateSets.get(index) == 0) {
            duplicateSets.set(index, duplicateSet);
            for (int i = 0; i < adjacencyList.get(index).size(); i++) {
                // Travese all members of the adjacency list to test if they belong to the current duplicate set
                mergeDuplicateSets(adjacencyList.get(index).get(i), duplicateSet, duplicateSets, adjacencyList);
            }
        }
    }

    private int getEditDistance(final String s1, final String s2) {
        if(s1 == null || s2 == null) {
            throw new PicardException("Attempt to compare two incomparable UMIs.  At least one of the UMIs was null.");
        }
        if(s1.length() != s2.length()) {
            throw new PicardException("Barcode " + s1 + " and " + s2 + " do not have matching lengths.");
        }
        int count = 0;
        for(int i = 0;i < s1.length();i++) {
            if(s1.charAt(i) != s2.charAt(i)) {
                count++;
            }
        }
        return count;
    }
}

