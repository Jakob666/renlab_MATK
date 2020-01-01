package SeqDataModel;

import java.util.Arrays;

/**
 * Size factor
 */
public class SizeFactor {
    private int[] ipGenesReads, inputGeneReads;

    /**
     * Constructor
     * @param ipGenesReads IP reads count for genes, shape 1 × individualNumber
     * @param inputGeneReads INPUT reads count for genes, shape 1 × individualNumber
     */
    public SizeFactor(int[] ipGenesReads, int[] inputGeneReads) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
    }

    private double calcSizeFactor(double readsCount1, double readsCount2) {
        readsCount2 = (readsCount2 <= 0.00001)? readsCount2 + 0.01: readsCount2;
        return readsCount1 / Math.sqrt(readsCount1 * readsCount2);
    }

    /**
     * IP data size factor of a single gene
     * @param ipRadsCount IP reads count
     * @param inputReadsCount INPUT reads count
     * @return IP size factor
     */
    private double calcIPSizeFactor(double ipRadsCount, double inputReadsCount) {
        return calcSizeFactor(ipRadsCount, inputReadsCount);
    }

    /**
     * INPUT data size factor of a single gene
     * @param inputReadsCount INPUT reads count
     * @param ipReadsCount IP reads count
     * @return INPUT size factor
     */
    private double calcINPUTSizeFactor(double inputReadsCount, double ipReadsCount) {
        return calcSizeFactor(inputReadsCount, ipReadsCount);
    }

    /**
     * global size factor
     * @param input true, if get INPUT data size factor; otherwise, IP data size factor
     * @return global size factor
     */
    public double getSizeFactor(boolean input) {
        double[] genesSizeFactor = new double[this.ipGenesReads.length];
        int ipReadsCount, inputReadsCount;
        for (int i=0; i<genesSizeFactor.length; i++) {
            ipReadsCount = this.ipGenesReads[i];
            inputReadsCount = this.inputGeneReads[i];
            if (input)
                genesSizeFactor[i] = calcINPUTSizeFactor(inputReadsCount, ipReadsCount);
            else
                genesSizeFactor[i] = calcIPSizeFactor(ipReadsCount, inputReadsCount);
        }

        Arrays.sort(genesSizeFactor);
        int mediaIndex = genesSizeFactor.length / 2;

        double sizeFactor;
        if (genesSizeFactor.length == 1) {
            return genesSizeFactor[0];
        }else if (genesSizeFactor.length % 2 == 0) {
            if (genesSizeFactor.length == 2)
                sizeFactor = (genesSizeFactor[0] + genesSizeFactor[1]) / 2;
            else
                sizeFactor = (genesSizeFactor[mediaIndex] + genesSizeFactor[mediaIndex+1]) / 2;
        } else {
            sizeFactor = genesSizeFactor[mediaIndex];
        }

        genesSizeFactor = null;
        return sizeFactor;
    }
}
