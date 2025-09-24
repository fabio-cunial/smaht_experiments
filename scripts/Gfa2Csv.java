import java.io.*;
import java.util.*;


public class Gfa2Csv {
    
    /**
     * Remark: the reference path is not used for computing counts.
     *
     * @param args
     */
    public static void main(String[] args) throws Exception {
        final String INPUT_GFA = args[0];
        final int MIN_COUNT = Integer.parseInt(args[1]);
        
        final String REF_PATH_NAME = "ref";
        
        int i;
        int count;
        String str, key;
        BufferedReader br;
        HashMap<String,Integer> counts;
        Map.Entry<String,Integer> entry;
        Iterator<Map.Entry<String,Integer>> iterator;
        String[] tokens, tokensPrime;
        
        
        // Computing node counts
        counts = new HashMap();
        br = new BufferedReader(new FileReader(INPUT_GFA));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)!='P') { str=br.readLine(); continue; }
            tokens=str.split("\t");
            if (tokens[1].equalsIgnoreCase(REF_PATH_NAME)) { str=br.readLine(); continue; }
            tokensPrime=tokens[2].split(",");
            for (i=0; i<tokensPrime.length; i++) {
                key=tokensPrime[i].substring(0,tokensPrime[i].length()-1);
                if (!counts.containsKey(key)) counts.put(key,Integer.valueOf(1));
                else counts.put(key,Integer.valueOf(counts.get(key).intValue()+1));
            }
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        System.out.println("Name,n_paths,Colour");
        iterator=counts.entrySet().iterator();
        while (iterator.hasNext()) {
            entry=iterator.next();
            count=entry.getValue().intValue();
            System.out.println(entry.getKey()+","+count+","+(count>=MIN_COUNT?"blue":"gray"));
        }
	}

}