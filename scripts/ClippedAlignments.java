import java.io.*;
import java.util.*;


/**
 * ----------> Also count reads with no clips but that have supplementary 
 * alignments to the same or different chrs, which we observed in ONT.
 */
public class ClippedAlignments {
    
    private static final int MIN_CLIP_LENGTH = 100;  // Arbitrary
    
    private static long nReads, nReadsWithClips;
    
    
    /**
     * The program reads from STDIN and writes to STDOUT and STDERR.
     *
     * @param args
     */
	public static void main(String[] args) throws Exception {
        final int CAPACITY = 100;  // Arbitrary

        boolean printRead;
        int i;
        int lastAlignment, leftClip, rightClip;
        long nAlignments;
        String str, readId, currentReadId, currentRead;
        BufferedReader stdin;
        BufferedWriter stdout, stderr;
        int[] tmpArray = new int[50000];
        String[] tokens;
        Alignment[] alignments;
        
        
        alignments = new Alignment[CAPACITY];
        for (i=0; i<alignments.length; i++) alignments[i] = new Alignment();
        stdin = new BufferedReader(new InputStreamReader(System.in));
        stdout = new BufferedWriter(new OutputStreamWriter(System.out));
        stderr = new BufferedWriter(new OutputStreamWriter(System.err));
        str=stdin.readLine(); currentReadId=""; currentRead=""; lastAlignment=-1; nAlignments=0; nReads=0; nReadsWithClips=0;
        while (str!=null) {
            nAlignments++;
            if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments. Total reads: "+nReads+" Reads with clipped alignments: "+nReadsWithClips+" ("+((100.0*nReadsWithClips)/nReads)+"%)");
            tokens=str.split("\t");
            leftClip=getSoftClip(tokens[5],true);
            rightClip=getSoftClip(tokens[5],false);
            readId=tokens[0];
            if (currentReadId.length()==0) {
                currentReadId=readId; currentRead=tokens[9];
                if (leftClip>=MIN_CLIP_LENGTH || rightClip>=MIN_CLIP_LENGTH) {
                    lastAlignment=0;
                    alignments[0].chr=tokens[2];
                    alignments[0].readStart=leftClip;
                    alignments[0].readEnd=tokens[9].length()-rightClip-1;
                    alignments[0].orientation=(Integer.parseInt(tokens[1])&16)!=0;
                }
                else lastAlignment=-1;
            }
            else if (!readId.equals(currentReadId)) {
                printRead=processAlignments(alignments,lastAlignment,currentRead.length(),stdout,tmpArray);
                if (printRead) {
                    stderr.write(">"+currentReadId); stderr.newLine();
                    stderr.write(currentRead); stderr.newLine();
                }
                // Next iteration
                currentReadId=readId; currentRead=tokens[9];
                if (leftClip>=MIN_CLIP_LENGTH || rightClip>=MIN_CLIP_LENGTH) {
                    lastAlignment=0;
                    alignments[0].chr=tokens[2];
                    alignments[0].readStart=leftClip;
                    alignments[0].readEnd=tokens[9].length()-rightClip-1;
                    alignments[0].orientation=(Integer.parseInt(tokens[1])&16)!=0;
                }
                else lastAlignment=-1;
            }
            else {
                if (leftClip>=MIN_CLIP_LENGTH || rightClip>=MIN_CLIP_LENGTH) {
                    lastAlignment++;
                    if (lastAlignment>=alignments.length) {
                        Alignment[] newArray = new Alignment[alignments.length<<1];
                        System.arraycopy(alignments,0,newArray,0,alignments.length);
                        for (i=alignments.length; i<newArray.length; i++) newArray[i] = new Alignment();
                        alignments=newArray;
                    }
                    alignments[lastAlignment].chr=tokens[2];
                    alignments[lastAlignment].readStart=leftClip;
                    alignments[lastAlignment].readEnd=tokens[9].length()-rightClip-1;
                    alignments[lastAlignment].orientation=(Integer.parseInt(tokens[1])&16)!=0;
                }
            }
            str=stdin.readLine();
        }
        printRead=processAlignments(alignments,lastAlignment,currentRead.length(),stdout,tmpArray);
        if (printRead) {
            stderr.write(">"+currentReadId); stderr.newLine();
            stderr.write(currentRead); stderr.newLine();
        }
        stdin.close(); stdout.close();
        System.err.println("Total reads: "+nReads);
        System.err.println("Reads with clipped alignments: "+nReadsWithClips+" ("+((100.0*nReadsWithClips)/nReads)+"%)");
	}
    
    
    /**
     * @param left check for a soft clip on the left (TRUE) or right (FALSE) 
     * side;
     * @return the length of the soft clip, or zero if no soft clip exists on
     * the specified side.
     */
    private static final int getSoftClip(String cigar, boolean left) {
        final int LENGTH = cigar.length();
        char c;
        int i;
        
        if (left) {
            for (i=0; i<LENGTH; i++) {
                c=cigar.charAt(i);
                if (!Character.isDigit(c)) return c=='S'?Integer.parseInt(cigar.substring(0,i)):0;
            }
        }
        else {
            if (cigar.charAt(LENGTH-1)!='S') return 0;
            for (i=LENGTH-2; i>=0; i--) {
                c=cigar.charAt(i);
                if (!Character.isDigit(c)) return Integer.parseInt(cigar.substring(i+1,LENGTH-1));
            }
        }
        return 0;
    }

    
    /**
     * @param tmpArray the procedure prints here a simple representation of 
     * `alignments` in read coordinates:
     *
     * -1: the read does not extend to this position (it is shorter);
     *  0: the read extends to this position, but no alignment covers the 
     *     position;
     *  1+X: a forward alignment to chrX covers this position (X>=1);
     * -1-X: a RC alignment to chrX covers this position (X>=1).
     *
     * @return TRUE iff there is a region >=500bp with no alignment at the 
     * beginning or end of the read.
     */
    private static final boolean processAlignments(Alignment[] alignments, int lastAlignment, int readLength, BufferedWriter bw, int[] tmpArray) throws IOException {
        final int QUANTUM = 10;  // Arbitrary
        final int MIN_UNMAPPED_LENGTH = 500;  // Arbitrary
        boolean found;
        int i, j;
        int numerator, denominator, minStart, maxEnd;
        
        nReads++;
        if (lastAlignment==-1) {
            // The read has no clipped alignment
            return false;
        }
        if (lastAlignment==0 && alignments[0].readEnd-alignments[0].readStart>readLength-MIN_CLIP_LENGTH) {
            System.err.println("ERROR: the read is almost fully covered by one clipped alignment whose clips are too short.");
            System.exit(1);
        }
        nReadsWithClips++;
        minStart=readLength; maxEnd=0;
        for (i=0; i<=lastAlignment; i++) {
            if (alignments[i].readStart<minStart) minStart=alignments[i].readStart;
            if (alignments[i].readEnd>maxEnd) maxEnd=alignments[i].readEnd;
        }
        if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
        Arrays.fill(tmpArray,-1);
        Arrays.fill(tmpArray,0,readLength>tmpArray.length?tmpArray.length:readLength,0);
        for (i=0; i<=lastAlignment; i++) alignments[i].print(tmpArray);
        for (i=0; i<tmpArray.length; i+=QUANTUM) {
            numerator=0; denominator=0;
            for (j=i; j<i+QUANTUM; j++) {
                if (j<tmpArray.length) { numerator+=tmpArray[j]; denominator++; }
            }
            if (denominator!=0) bw.write((numerator/denominator)+",");
        }
        bw.newLine();
        return minStart>=MIN_UNMAPPED_LENGTH || maxEnd<readLength-MIN_UNMAPPED_LENGTH;
    }
    
    
    /**
     * @return 1+X where X>=1 is an integer assigned to `chr`.
     */
    private static int chr2int(String chr) {
        final int NONCANONICAL_CHR = 30;  // Arbitrary
        
        if (chr.equalsIgnoreCase("chrX")) return 1+23;
        else if (chr.equalsIgnoreCase("chrY")) return 1+24;
        else if (chr.equalsIgnoreCase("chrM")) return 1+25;
        else return chr.length()>5?NONCANONICAL_CHR:(1+Integer.parseInt(chr.substring(3)));
    }
    
    
    private static class Alignment implements Comparable {
        public boolean orientation;
        public String chr;
        public int readStart, readEnd;  // Zero-based
        
        public int compareTo(Object other) {
            Alignment otherAlignment = (Alignment)other;
            if (readStart<otherAlignment.readStart) return -1;
            else if (readStart>otherAlignment.readStart) return 1;
            if (readEnd<otherAlignment.readEnd) return -1;
            else if (readEnd>otherAlignment.readEnd) return 1;
            return 0;
        }
        
        /**
         * Remark: the procedure overwrites existing values in `out`, if any.
         */
        public void print(int[] out) {
            final int MAX_POS = Math.min(readEnd,out.length-1);
            final int VALUE = orientation?chr2int(chr):-chr2int(chr);
            
            for (int i=readStart; i<=MAX_POS; i++) out[i]=VALUE;
        }
    }

}