package SeqDataModel;

import java.util.Arrays;

/**
 *                                                                              X_ip,i
 * Size factor of one sample. Suppose n genes, size-factor_ip = median --------------------------
 *                                                                      (X_ip,i * X_input,i)^0.5
 *
 *                                                                             X_input,i
 *                                             size_factor_input = median --------------------------
 *                                                                         (X_ip,i * X_input,i)^0.5
 *                            where i from 1 to n
 */
public class SizeFactor {
    private int geneNumber;
    private int[] ipGenesReads, inputGeneReads;

    /**
     * Constructor
     * @param ipGenesReads IP reads count for genes, shape 1 × geneNumber
     * @param inputGeneReads INPUT reads count for genes, shape 1 × geneNumber
     */
    public SizeFactor(int[] ipGenesReads, int[] inputGeneReads) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.geneNumber = ipGenesReads.length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
    }

    private double calcSizeFactor(double readsCount1, double readsCount2) {
        readsCount1 = (readsCount1 <= 0.00001)? readsCount1 + 0.01: readsCount1;
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
        double[] genesSizeFactor = new double[this.geneNumber];
        int ipReadsCount, inputReadsCount;
        // size factor of all genes
        for (int i=0; i<this.geneNumber; i++) {
            ipReadsCount = this.ipGenesReads[i];
            inputReadsCount = this.inputGeneReads[i];
            if (input)
                genesSizeFactor[i] = calcINPUTSizeFactor(inputReadsCount, ipReadsCount);
            else
                genesSizeFactor[i] = calcIPSizeFactor(ipReadsCount, inputReadsCount);
        }

        Arrays.sort(genesSizeFactor);
        int mediaIndex = genesSizeFactor.length / 2;

        // use median as size factor of this individual
        double sizeFactor;
        if (this.geneNumber == 1)
            return genesSizeFactor[0];
        if (this.geneNumber % 2 == 0)
            sizeFactor = genesSizeFactor[mediaIndex+1];
        else
            sizeFactor = (genesSizeFactor[mediaIndex-1] + genesSizeFactor[mediaIndex]) / 2;
        genesSizeFactor = null;
        return sizeFactor;
    }
}
