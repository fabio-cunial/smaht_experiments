import java.io.*;
import java.util.*;


public class ClippedAlignments {
    
    private static long nReads, nReadsWithClips;
    
    
    /**
     * @param args
     */
	public static void main(String[] args) throws Exception {
		final String INPUT_SAM = args[0];
        final String OUTPUT_CSV = args[1];
        
        final int CAPACITY = 100;  // Arbitrary
        final int MIN_CLIP_LENGTH = 100;  // Arbitrary

        int i;
        int lastAlignment, leftClip, rightClip;
        long nAlignments;
        String str, readId, currentReadId;
        BufferedReader br;
        BufferedWriter bw;
        int[] tmpArray = new int[30000];
        String[] tokens;
        Alignment[] alignments;
        
        
        alignments = new Alignment[CAPACITY];
        for (i=0; i<alignments.length; i++) alignments[i] = new Alignment();
        bw = new BufferedWriter(new FileWriter(OUTPUT_CSV));
        br = new BufferedReader(new FileReader(INPUT_SAM));
        str=br.readLine(); currentReadId=""; lastAlignment=-1; nAlignments=0; nReads=0; nReadsWithClips=0;
        while (str!=null) {
            nAlignments++;
            if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments...");
            tokens=str.split("\t");
            leftClip=getSoftClip(tokens[5],true);
            rightClip=getSoftClip(tokens[5],false);
            readId=tokens[0];
            if (currentReadId.length()==0) {
                currentReadId=readId;
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
                processAlignments(alignments,lastAlignment,bw,tmpArray);
                // Next iteration
                currentReadId=readId;
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
            str=br.readLine();
        }
        processAlignments(alignments,lastAlignment,bw,tmpArray);
        br.close(); bw.close();
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
     * Prints a simple representation of `alignments` in read coordinates.
     */
    private static final void processAlignments(Alignment[] alignments, int lastAlignment, BufferedWriter bw, int[] tmpArray) throws IOException {
        final int QUANTUM = 10;  // Arbitrary
        int i, j;
        int numerator, denominator;
        
        nReads++;
        if (lastAlignment>=0) nReadsWithClips++;
        if (lastAlignment>0) Arrays.sort(alignments,0,lastAlignment+1);
        Arrays.fill(tmpArray,0);
        for (i=0; i<=lastAlignment; i++) alignments[i].print(tmpArray);
        for (i=0; i<tmpArray.length; i+=QUANTUM) {
            numerator=0; denominator=0;
            for (j=i; j<i+QUANTUM; j++) {
                if (j<tmpArray.length) { numerator+=tmpArray[j]; denominator++; }
            }
            if (denominator!=0) bw.write((numerator/denominator)+",");
        }
        bw.newLine();
    }
    
    
    private static int chr2int(String chr) {    
        if (chr.equalsIgnoreCase("chrX")) return 23;
        else if (chr.equalsIgnoreCase("chrY")) return 24;
        else if (chr.equalsIgnoreCase("chrM")) return 25;
        else return chr.length()>5?Integer.MAX_VALUE:Integer.parseInt(chr.substring(3));
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