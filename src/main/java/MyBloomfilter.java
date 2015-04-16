import ErrorCorrection.Kmer_general;
import com.skjegstad.utils.BloomFilter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static com.skjegstad.utils.BloomfilterBenchmark.printStat;

public class MyBloomfilter {
/**This class is to input the kmers from different fasta files and add the kmers of the
 *  file into bloomfilter. if the kmer already in the bF database, then ? */
    /**
     Parameters:
     args[0]:is the input dir for fasta files including the initial fas files and query fas files, includes 2 sub folders
     args[1]: is kmer length here is 26
     args[2]:is the output information for the kmers including the fasta file kmers and the true false kmers not
     included in the bloom filter.
     inputDir: is the directory which hold the input fasta files
     singleFasta: is the fas file in the input directory.
     outputDir: is the directory which hold the output log files
     outputDir2: is the directory which hold the true false kmer of the fasta file
     fileKmersAll:
     */
    public static void main (String[] args) throws IOException {
        // inputDir includs two directories,one is the ref fastas, the other is the query fastas
        File inputDir = new File(args[0]);
        //System.out.println("args[0]: " + args[0]);
        if (args.length != 4){
            System.out.println("usage message: java -jar inputDir kmerSize outputDir" + "\n" +
                    "java required version is 1.8");
        }
        File[] in_files = inputDir.listFiles();
        //String maxV = "" + args[4];
        BloomFilter<String> bf = new BloomFilter<String>(Double.parseDouble(args[3]), Integer.MAX_VALUE);
        System.out.println("k is optimal number of hash functions: " + bf.getK());
        int kmerCount = 0;
        int reffileCounter = 0;
        int queryfileCounter = 0;
        long start_add = System.currentTimeMillis();
        Map<String, Integer> insideKmerInfor_ini = new HashMap<String, Integer>();
        Map<String, Integer> insideKmerInfor_que = new HashMap<String, Integer>();
        for (int i = 0; i < in_files.length; i++){
            /** input directory have two sub-directories, one have all the initial(ref) fasta fils, the other one have all
             * the query fasta files*/
            //deal with initial fasta files
            if (in_files[i].getName().equals("initial_database")) {
                File[] in_fastas = in_files[i].listFiles(); //in_fastas is the folder which have all initial fastas
                String[] singleFasta = new String[in_fastas.length];
                String[] initial_outputlogfile = new String[in_fastas.length];
                String[] initial_outputfalselogfile = new String[in_fastas.length];
                List<Kmer_general> singleFileKmer_ini;
                String insideKmerLog_ini = args[2] + "/" + "databaseFas_insideKmer.log";
                for (int nu_fastas = 0; nu_fastas < in_fastas.length; nu_fastas++) {
                    List<Kmer_general> falseKmer_ini = new ArrayList<Kmer_general>();
                    List<Kmer_general> insideKmer_ini = new ArrayList<Kmer_general>();
                    reffileCounter = reffileCounter + 1;
                    System.out.println("\n" + "Ref File counter ************ " + reffileCounter + "*************" + "\n");
                    singleFasta[nu_fastas] = in_fastas[nu_fastas].getName();
                    //String outputDir = args[2] + "logfile" + Integer.toString(i)+".log" ;
                    File output_initial = new File(args[2] + "/"+ "initial_output");
                    File outputlog_initial = new File(output_initial.getPath() + "/" + "AllKmerLog");
                    File outputfalselog_initial = new File(output_initial.getPath() + "/" + "falseKmerLog");
                    outputlog_initial.mkdirs();
                    outputfalselog_initial.mkdirs();
                    initial_outputlogfile[nu_fastas] = outputlog_initial.getPath()+ "/" + singleFasta[nu_fastas] + ".log";
                    initial_outputfalselogfile[nu_fastas] = outputfalselog_initial.getPath()+ "/" + singleFasta[nu_fastas] + ".false.log";
                    System.out.println("mybloomfilter.java : the input file is " + singleFasta[nu_fastas]);
                    ReDataset ds_ref = new ReDataset(in_fastas[nu_fastas].getPath());
                    /**setK method is public in the original author's code*/
                    ds_ref.setK(Integer.parseInt(args[1]));
                    ds_ref.calculateKMersAndKCounts();
                    singleFileKmer_ini = ds_ref.getFileKmers();
                    kmerCount = kmerCount + ds_ref.getFileKmers().size();
                    ds_ref.printout(initial_outputlogfile[nu_fastas]);
                    /**-----------------------------------------------------------*/
                    /**
                     * Returns the value chosen for K.<br />
                     * <br />
                     * K is the optimal number of hash functions based on the size
                     * of the Bloom filter and the expected number of inserted elements.
                     *
                     * @return optimal k.
                     */
                    /**add elements. 1. check for existing elements with contains()
                     * if the kmer does not in, then add it into BF, else ?*/
                    int falseCounter = 0;
                    int insideKmerCounter = 0;
                    for (int j = 0; j < singleFileKmer_ini.size(); j++) {
                        boolean in = bf.contains(singleFileKmer_ini.get(j).getNucl());
                        if (in != true) {
                            bf.add(ds_ref.getFileKmers().get(j).getNucl());
                            falseCounter = falseCounter + 1;
                            falseKmer_ini.add(singleFileKmer_ini.get(j));
                        }else if (in == true){
                            insideKmerCounter = insideKmerCounter + 1;
                            insideKmer_ini.add(singleFileKmer_ini.get(j));
                        }
                    }
                    insideKmerInfor_ini.put(singleFasta[nu_fastas],insideKmerCounter);
                    System.out.println("not inside kmer for ref file: " + falseCounter);
                    System.out.println("inside kmer for ref file: " + insideKmerCounter);
                    Tools.printoutLog(initial_outputfalselogfile[nu_fastas],falseKmer_ini);
                }
                Tools.printoutInsideKmerLog(insideKmerLog_ini, insideKmerInfor_ini);
             //deal with query fasta files
            } else if (in_files[i].getName().equals("query")){
                long start_queryadd = System.currentTimeMillis();
                File[] query_fastas = in_files[i].listFiles(); //query_fastas is the folder which have all query fasta files
                String[] singleFasta_q = new String[query_fastas.length];
                String[] query_outputlogfile = new String[query_fastas.length];
                String[] query_outputfalselogfile = new String[query_fastas.length];
                List<Kmer_general> singleFileKmer_q;
                String insideKmerLog_query = args[2] + "/" + "queryFas_insideKmer.log";
                for (int nu_qfastas = 0; nu_qfastas < query_fastas.length; nu_qfastas++) {
                    List<Kmer_general> falseKmer_query = new ArrayList<Kmer_general>();
                    List<Kmer_general> insideKmer_query = new ArrayList<Kmer_general>();
                    queryfileCounter = queryfileCounter + 1;
                    System.out.println("\n" + "Query File counter ************ " + queryfileCounter + "*************" + "\n");
                    singleFasta_q[nu_qfastas] = query_fastas[nu_qfastas].getName();
                    File output_query = new File(args[2] + "/" + "query_output");
                    File outputlog_query = new File(output_query + "/" + "AllKmerLog");
                    File outputfalselog_query = new File(output_query + "/" + "falseKmerLog");
                    outputlog_query.mkdirs();
                    outputfalselog_query.mkdirs();
                    query_outputlogfile[nu_qfastas] = outputlog_query.getPath() + "/"+ singleFasta_q[nu_qfastas] + ".log";
                    query_outputfalselogfile[nu_qfastas] = outputfalselog_query.getPath()+ "/" + singleFasta_q[nu_qfastas] + ".false.log";
                    System.out.println("mybloomfilter.java : the input file is " + singleFasta_q[nu_qfastas]);
                    ReDataset ds_ref = new ReDataset(query_fastas[nu_qfastas].getPath());
                    /**setK method is public in the original author's code*/
                    ds_ref.setK(Integer.parseInt(args[1]));
                    ds_ref.calculateKMersAndKCounts();
                    singleFileKmer_q = ds_ref.getFileKmers();
                    kmerCount = kmerCount + ds_ref.getFileKmers().size();
                    ds_ref.printout(query_outputlogfile[nu_qfastas]);
                    /**-----------------------------------------------------------*/
                    /**
                     * Returns the value chosen for K.<br />
                     * <br />
                     * K is the optimal number of hash functions based on the size
                     * of the Bloom filter and the expected number of inserted elements.
                     *
                     * @return optimal k.
                     */
                    /**add elements. 1. check for existing elements with contains()
                     * if the kmer does not in, then add it into BF, else ?*/
                    int falseCounter = 0;
                    int insideKmerCounter_q = 0;
                    for (int j = 0; j < singleFileKmer_q.size(); j++) {
                        boolean in = bf.contains(singleFileKmer_q.get(j).getNucl());
                        if (in != true) {
                            bf.add(ds_ref.getFileKmers().get(j).getNucl());
                            falseCounter = falseCounter + 1;
                            falseKmer_query.add(singleFileKmer_q.get(j));
                        }else if (in == true){
                            insideKmerCounter_q = insideKmerCounter_q + 1;
                            insideKmer_query.add(singleFileKmer_q.get(j));
                        }
                    }
                    insideKmerInfor_que.put(singleFasta_q[nu_qfastas], insideKmerCounter_q);
                    System.out.println("Not inside kmer for query file: " + falseCounter);
                    System.out.println("Inside kmer for query file: " + insideKmerCounter_q);
                    Tools.printoutLog(query_outputfalselogfile[nu_qfastas], falseKmer_query);
                }
                Tools.printoutInsideKmerLog(insideKmerLog_query,insideKmerInfor_que);
                long end_queryadd = System.currentTimeMillis();
                System.out.print("The time used in add query kmers into bloom filter is: ");
                printStat(start_queryadd, end_queryadd);
            }
        }
            /**-------------------------------------------------------------*/
        long end_add = System.currentTimeMillis();
        System.out.print("The time used in add all kmers into bloom filter is: ");
        printStat(start_add, end_add);
    }
}
