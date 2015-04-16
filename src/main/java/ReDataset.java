import ErrorCorrection.DataSet;
import ErrorCorrection.Kmer_general;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/*
 * Created by yqb7 on 3/17/2015.
 */

public class ReDataset extends DataSet{
    static ArrayList<Kmer_general> myKmers;
    public ReDataset( String inputfile) throws IOException {
        super(inputfile);
    }

    /**This method is to get the kmers for single fas file. variable myKmers
    includes kmer string and kmers counts*/
    public ArrayList<Kmer_general> getFileKmers( ) {
        /**AllKmers in the Pavel's code means single file's kmers*/
        myKmers = super.getAllKmers();
        return myKmers;
    }

    /**This method is to output the */
    public void printout(String fileKmer) throws IOException {
        String outputfile = fileKmer;
        FileWriter outputStream = null;
        try {
            outputStream = new FileWriter(new File(outputfile));
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < myKmers.size(); i++ ){
            outputStream.write(myKmers.get(i).getNucl() + " " +
                    Integer.toString(myKmers.get(i).kCount) + "\n");
        }
        outputStream.close();
    }
}

