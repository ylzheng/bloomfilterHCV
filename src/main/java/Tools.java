import ErrorCorrection.Kmer_general;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Created by yqb7 on 3/30/2015.
 */
public class Tools {

     /** This method is used to print the false kmer log file*/
    public static void printoutLog(String bfLog, List<Kmer_general> false_Kmer) throws IOException {
        FileWriter outputStream = null;
        try {
            outputStream = new FileWriter(new File(bfLog));
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < false_Kmer.size(); i++ ){
            outputStream.write(false_Kmer.get(i).getNucl() + " " +
                    Integer.toString(false_Kmer.get(i).kCount) +"\n");
        }
        outputStream.close();
    }

    public static void printoutInsideKmerLog(String bfInKmerLog, Map<String, Integer> inKmerlog_q) throws IOException {
        FileWriter outputStream = null;
        try {
            outputStream = new FileWriter(new File(bfInKmerLog));
        } catch (IOException e) {
            e.printStackTrace();
        }
        outputStream.write("FileName         #Kmers in BloomFilter" + "\n");
        Iterator<?> it = inKmerlog_q.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry<String, Integer> pairs = (Map.Entry<String, Integer>)it.next();
            outputStream.write(pairs.getKey() + "    " + pairs.getValue() + "\n");
        }
        outputStream.close();
    }
}
