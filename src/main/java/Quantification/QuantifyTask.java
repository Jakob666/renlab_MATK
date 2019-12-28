package Quantification;

public interface QuantifyTask {
    Runnable getTask(int[][] sampleIPReads, int[][] sampleINPUTReads, double[][] sampleIPExpectation, double[][] sampleINPUTExpectation, int geneIdx);
}
