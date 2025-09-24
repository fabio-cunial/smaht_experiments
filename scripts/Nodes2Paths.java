import java.io.*;
import java.util.*;


public class Nodes2Paths {
    
    /**
     * @param args 
     * 0: input GFA;
     * 1..: the node IDs through which paths should pass.
     */
    public static void main(String[] args) throws Exception {
        final String INPUT_GFA = args[0];
        
        int i, j;
        int nodeID;
        String str;
        StringBuilder sb;
        BufferedReader br;
        int[] queryNodes;
        String[] tokens;
        
        sb = new StringBuilder();
        queryNodes = new int[args.length-1];
        for (i=1; i<args.length; i++) queryNodes[i-1]=Integer.parseInt(args[i]);
        br = new BufferedReader(new FileReader(INPUT_GFA));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)!='P') { str=br.readLine(); continue; }
            tokens=str.split("\t");
            tokens=tokens[2].split(",");
            sb.delete(0,sb.length());
            for (i=0; i<tokens.length; i++) {
                nodeID=Integer.parseInt(tokens[i].substring(0,tokens[i].length()-1));
                for (j=0; j<queryNodes.length; j++) {
                    if (queryNodes[j]==nodeID) sb.append(tokens[i]+",");
                }
            }
            if (sb.length()!=0) System.out.println(sb.toString());
            str=br.readLine();
        }
        br.close();
	}

}