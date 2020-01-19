package Quantification;

import SeqDataModel.BackgroundExpression;

import java.util.ArrayList;

/**
 * Quantify the methylation level with component-wise MH sampling
 * reference https://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
 */
public class MethylationQuantification {
    private QuantificationModel quantificationModel;
    private int peakNumber, samplingTime, burnIn, threadNumber;
    private int[][] ipReads, inputReads;
    private double a, b, c, d, w, k, k1, k2;
    private String tmpDir, outputFile;

    public MethylationQuantification(int[][] ipReads, int[][] inputReads,
                                     double ipShape, double ipScale, double inputShape, double inputScale,
                                     double methLevelShape1, double methLevelShape2,
                                     double expansionEffect1, double expansionEffect2,
                                     int samplingTime, int burnIn, int threadNumber, String tmpDir, String outputFile) {
        assert samplingTime > burnIn;
        // assert arrays contains same sample and peak number
        assert inputReads.length == ipReads.length && inputReads[0].length == ipReads[0].length;
        this.peakNumber = inputReads[0].length;
        this.inputReads = inputReads;
        this.ipReads = ipReads;
        this.a = ipShape;
        this.b = ipScale;
        this.c = inputShape;
        this.d = inputScale;
        this.w = methLevelShape1;
        this.k = methLevelShape2;
        this.k1 = expansionEffect1;
        this.k2 = expansionEffect2;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;
        this.tmpDir = tmpDir;
        this.outputFile = outputFile;
    }

    /**
     * API function
     */
    public void runExecution() {
        this.initialize();
    }

    /**
     * get quantification result
     * @return quantification result
     */
    public ArrayList<String> getQuantifyResult() {
        ArrayList<String> result = new ArrayList<>();
        this.quantificationModel.getQuantifyResult().entrySet().stream().sorted((o1, o2) -> o1.getKey() - o2.getKey()).forEach(x -> result.add(x.getValue()));
        return result;
    }

    /**
     * initialize samplers with given parameters and quantification model
     */
    private void initialize() {
        MethylationLevelSampler methLevelSampler = new MethylationLevelSampler(this.w, this.k);
        OverdispersionSampler inputOverdispersionSampler = new OverdispersionSampler(this.c, this.d);
        OverdispersionSampler ipOverdispersionSampler = new OverdispersionSampler(this.a, this.b);
        ExpansionEffectSampler expansionEffectSampler = new ExpansionEffectSampler(this.k1, this.k2);
        BackgroundExpression be = new BackgroundExpression(this.ipReads, this.inputReads);
        // background expression expectation and standard deviation, shape 1 × peakNumber
        double[] backgroundExpressionMean = be.geneBackgroundExp();
        double[] backgroundExpressionStd = be.geneExpressionStd();

        // for each peak, initialize a background expression sampler. cancel background expression sampling
        BackgroundExpressionSampler[] geneBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.peakNumber];
        BackgroundExpressionSampler sampler;
        for (int i=0; i<this.peakNumber; i++) {
            // shape and scale parameters of log-normal distribution
            double scale = Math.log(backgroundExpressionMean[i]);
            double shape = backgroundExpressionStd[i];
            sampler = new BackgroundExpressionSampler(scale, shape);
            // random init background expression from corresponding distribution
            geneBackgroundExpressionSamplers[i] = sampler;
        }

        SamplingRecords samplingRecords = new SamplingRecords();
        double[] methylationLevels = new double[this.peakNumber];
        double[] expansionEffectValues = new double[this.samplingTime];
        double[] ipOverdispersions = new double[this.samplingTime];
        double[] inputOverdispersions = new double[this.samplingTime];
        expansionEffectValues[0] = 1;
        ipOverdispersions[0] = 0.1;
        inputOverdispersions[0] = 0.1;
        for (int i=0; i<this.peakNumber; i++) {
            methylationLevels[i] = 0.1;
        }
        samplingRecords.setMethylationLevels(methylationLevels);
        samplingRecords.setBackgroundExpressions(backgroundExpressionMean);
        samplingRecords.setExpansionEffects(expansionEffectValues);
        samplingRecords.setIpOverdispersions(ipOverdispersions);
        samplingRecords.setInputOverdispersions(inputOverdispersions);

        this.quantificationModel = new QuantificationModel(this.ipReads, this.inputReads, methLevelSampler, expansionEffectSampler,
                                                           ipOverdispersionSampler, inputOverdispersionSampler, geneBackgroundExpressionSamplers,
                                                           samplingRecords, this.threadNumber, this.samplingTime, this.burnIn,
                                                           this.tmpDir, this.outputFile);
    }

    /**
     * run execution and record into file (if output file is not null)
     */
    public void runModel() {
        this.quantificationModel.calcSampleSizeFactors();
        this.quantificationModel.iteration();
        this.quantificationModel.quantify();
        if (this.outputFile!=null)
            this.quantificationModel.output();
    }
}
