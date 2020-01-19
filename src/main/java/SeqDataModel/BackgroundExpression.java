package SeqDataModel;

/**
 * Background expression in the i th peak region, suppose there are m samples
 *               1        X_input,i,j
 *      exp_i = --- sum(---------------)    j from 1 to m
 *               m       size-factor_j
 *
 * where X_input,i,j is the reads count in i th peak region of the j th sample
 *       size-factor_j is the j th sample's size factor
 */
public class BackgroundExpression {
    private int individualNumber, geneNumber;
    private int[][] ipGenesReads, inputGeneReads;
    private double[] inputSizeFactors = null;

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape individual number × gene number
     * @param inputGeneReads INPUT gene reads, shape individual number × gene number
     */
    public BackgroundExpression(int[][] ipGenesReads, int[][] inputGeneReads) {
        // has same individual number
        assert ipGenesReads.length == inputGeneReads.length;
        this.individualNumber = ipGenesReads.length;
        this.geneNumber = ipGenesReads[0].length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
    }

    /**
     * calculate global IP size factor for each individual
     * @return IP size factor, shape 1 × individualNumber
     */
    public double[] getGlobalIPSizeFactor() {
        double[] individualSizeFactors = new double[this.individualNumber];  // shape 1 × individualNumber
        for (int i=0; i<this.individualNumber; i++) {
            // reads count of genes in IP and INPUT data of individual i
            int[] individualIPReads = this.ipGenesReads[i], individualINPUTReads = this.inputGeneReads[i];
            SizeFactor sf = new SizeFactor(individualIPReads, individualINPUTReads);
            double sizeFactor = sf.getSizeFactor(false);   // individual i's size factor
            sf = null;
            individualSizeFactors[i] = sizeFactor;
        }

        return individualSizeFactors;
    }

    /**
     * calculate global INPUT size factor for each individual
     * @return IP size factor, shape 1 × individualNumber
     */
    public double[] getGlobalINPUTSizeFactor() {
        double[] individualSizeFactors = new double[this.individualNumber];  // shape 1 × individualNumber
        for (int i=0; i<this.individualNumber; i++) {
            // reads count of genes in IP and INPUT data of individual i
            int[] individualIPReads = this.ipGenesReads[i], individualINPUTReads = this.inputGeneReads[i];
            SizeFactor sf = new SizeFactor(individualIPReads, individualINPUTReads);
            double sizeFactor = sf.getSizeFactor(true);   // individual i's size factor
            sf = null;
            individualSizeFactors[i] = sizeFactor;
        }

        return individualSizeFactors;
    }

    /**
     * calculate gene expression for each gene
     * @return gene background expression, shape 1 × geneNumber
     */
    public double[] geneBackgroundExp() {
        this.inputSizeFactors = this.getGlobalINPUTSizeFactor();   // shape 1 × individualNumber
        double[] genesBkgExp = new double[this.geneNumber];  // shape 1 × geneNumber

        for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
            double individualSizeFactor = this.inputSizeFactors[individualIdx];
            int[] individualGeneReads = this.inputGeneReads[individualIdx];  // one sample INPUT data observed reads

            for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
                genesBkgExp[geneIdx] += individualGeneReads[geneIdx] / individualSizeFactor;
            }
        }

        for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
            genesBkgExp[geneIdx] = genesBkgExp[geneIdx] / this.individualNumber;
        }

        return genesBkgExp;
        // only for simulation data test
//        double[] res = new double[this.geneNumber];
//        for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
//            res[geneIdx] = 500;
//        }
//        return res;
    }

    /**
     * calculate gene expression for each gene
     * @return gene background expression, shape 1 × geneNumber
     */
    public double[] geneBackgroundExp(double[] individualINPUTSizeFactors) {
        // shape 1 × individualNumber
        this.inputSizeFactors = individualINPUTSizeFactors;
        // shape 1 × geneNumber
        double[] genesBkgExp = new double[this.geneNumber];

        for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
            double individualSizeFactor = this.inputSizeFactors[individualIdx];
            int[] individualGeneReads = this.ipGenesReads[individualIdx];

            for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
                genesBkgExp[geneIdx] += individualGeneReads[geneIdx] / individualSizeFactor;
            }
        }

        for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
            genesBkgExp[geneIdx] = genesBkgExp[geneIdx] / this.individualNumber;
        }

        return genesBkgExp;
        // only for simulation data test
//        double[] res = new double[this.geneNumber];
//        for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
//            res[geneIdx] = 500;
//        }
//        return res;
    }

    /**
     * gene background expression standard  deviation
     * @return standard deviation, shape 1 × geneNumber
     */
    public double[] geneExpressionStd() {
         double[] res = new double[this.geneNumber];
         for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
             res[geneIdx] = 0.1;
         }
         return res;
    }
}
