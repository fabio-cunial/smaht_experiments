import java.io.*;
import java.util.*;
import java.util.zip.*;


public class AnalyzeClippedAlignments {
    
    /**
     * @param args
     */
	public static void main(String[] args) throws Exception {
        final String INPUT_CSV_GZ = args[0];
        final String OUTPUT_PREFIX = args[1];
        
        final int NONCANONICAL_CHR = 30;
        
        int i, j;
        int last, value;
        long nReads;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        int[] values;
        long[] chromosome_counts;
        long[][] pair_counts;
        
        // Counting
        chromosome_counts = new long[31];
        pair_counts = new long[31][31];
        values = new int[3001];
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_CSV_GZ))));
        str=br.readLine(); nReads=0;
        while (str!=null) {
            tokens=str.split(",");
            for (i=0; i<tokens.length; i++) {
                value=Integer.parseInt(tokens[i]);
                if (value<0) value=-value;
                values[i]=value;
            }
            Arrays.sort(values);
            last=0;
            for (i=1; i<values.length; i++) {
                if (values[i]==values[last]) continue;
                values[++last]=values[i];
            }
            for (i=0; i<=last; i++) {
                if (values[i]<=25 || values[i]==NONCANONICAL_CHR) chromosome_counts[values[i]]++;
            }
            for (i=0; i<=last; i++) {
                for (j=i+1; j<=last; j++) {
                    if ((values[i]<=25 || values[i]==NONCANONICAL_CHR) && (values[j]<=25 || values[j]==NONCANONICAL_CHR)) pair_counts[values[i]][values[j]]++;
                }
            }
            nReads++;
            if (nReads%10000==0) System.err.println("Analyzed "+nReads+" reads...");
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_chromosome_counts.csv"));
        for (i=0; i<chromosome_counts.length; i++) bw.write(i+","+chromosome_counts[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_pair_counts.csv"));
        for (i=0; i<pair_counts.length; i++) {
            for (j=0; j<pair_counts[i].length; j++) bw.write(pair_counts[i][j]+",");
            bw.newLine();
        }
        bw.close();
	}

}