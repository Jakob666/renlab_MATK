package SeqDataModel;

import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.util.Arrays;

public class BackgroundExpression {
    private int[][] ipGenesReads, inputGeneReads;

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param inputGeneReads INPUT gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     */
    public BackgroundExpression(int[][] ipGenesReads, int[][] inputGeneReads) {
        // has same individual number
        assert ipGenesReads.length == inputGeneReads.length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
    }

    /**
     * calculate global IP size factor for each individual
     * @return IP size factor, shape 1 × individualNumber
     */
    private double[] getGlobalSizeFactor() {
        // shape 1 × individualNumber
        double[] individualSizeFactors = new double[ipGenesReads.length];
        for (int i=0; i<individualSizeFactors.length; i++) {
            // reads count of genes in IP and INPUT data of individual i
            int[] individualIPReads = this.ipGenesReads[i], individualINPUTReads = this.inputGeneReads[i];
            SizeFactor sf = new SizeFactor(individualIPReads, individualINPUTReads);
            // individual i's size factor
            double sizeFactor = sf.getSizeFactor(false);
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
        // shape 1 × individualNumber
        double[] individualSizeFactors = this.getGlobalSizeFactor();
        // shape 1 × geneNumber
        double[] genesBkgExp = new double[this.ipGenesReads[0].length];

        int individualNumber = individualSizeFactors.length;
        for (int individualIdx=0; individualIdx<individualNumber; individualIdx++) {
            double individualSizeFactor = individualSizeFactors[individualIdx];
            int[] individualGeneReads = this.ipGenesReads[individualIdx];

            for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
                genesBkgExp[geneIdx] += individualGeneReads[geneIdx] / individualSizeFactor;
            }
        }

        for (int geneIdx=0; geneIdx<genesBkgExp.length; geneIdx++) {
            genesBkgExp[geneIdx] = genesBkgExp[geneIdx] / individualNumber;
        }

        return genesBkgExp;
        // TODO:当作用于实际数据时改用上一行被注释的代码
//        double[] res = new double[genesBkgExp.length];
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
        int individualNumber = this.ipGenesReads.length;
        int geneNumber = this.ipGenesReads[0].length;
//        if (individualNumber == 1)
//            return new double[geneNumber];

        // genes' background expression expectation, shape 1 × geneNumber
//        double[] genesBkgExpMean = this.geneBackgroundExp();
        // genes' background expression standard deviation, shape 1 × geneNumber
//        double[] genesBkgExpStd = new double[this.ipGenesReads[0].length];

//        double distance;
//        for (int[] individualGeneIPReads: this.ipGenesReads) {
//            for (int geneIdx=0; geneIdx<geneNumber; geneIdx++) {
//                distance = individualGeneIPReads[geneIdx] - genesBkgExpMean[geneIdx];
//                genesBkgExpStd[geneIdx] += Math.pow(distance, 2);
//            }
//        }

//        for (int geneIdx=0; geneIdx<genesBkgExpStd.length; geneIdx++) {
//            genesBkgExpStd[geneIdx] = Math.sqrt(genesBkgExpStd[geneIdx] / individualNumber);
//        }
//        genesBkgExpMean = null;

//        return genesBkgExpStd;

        // TODO:当作用于实际数据时改用上方的代码
        double[] res = new double[this.ipGenesReads[0].length];
        for (int geneIdx=0; geneIdx<this.ipGenesReads[0].length; geneIdx++) {
            res[geneIdx] = 0.1;
        }
        return res;
    }
}
