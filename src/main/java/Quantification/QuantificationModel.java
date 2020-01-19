package Quantification;

import ProbabilityCalculation.ProbabilityCalculator;
import SeqDataModel.SizeFactor;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.ReentrantLock;


public class QuantificationModel {
    private int individualNumber, peakNumber, threadNumber, samplingTime, burnIn;
    private int[][] ipReads, inputReads;    // shape individualNumber × peakNumber
    private double ipOverdispersion, inputOverdispersion, expansionEffect;
    private String tmpDir, outputFile;
    private HashMap<Integer, String> quantifyResult;
    private MethylationLevelSampler methylationLevelSampler;
    private ExpansionEffectSampler expansionEffectSampler;
    private OverdispersionSampler ipOverdispersionSampler, inputOverdispersionSampler;
    private BackgroundExpressionSampler[] geneBackgroundExpressionSamplers;
    private SamplingRecords samplingRecords;
    private ExecutorService executorService;
    private ReentrantLock lock;
    private DecimalFormat df = new DecimalFormat("0.0000");

    public QuantificationModel(int[][] ipReads, int[][] inputReads,
                               MethylationLevelSampler methLevelSampler, ExpansionEffectSampler expansionEffectSampler,
                               OverdispersionSampler ipOverdispersionSampler, OverdispersionSampler inputOverdispersionSampler,
                               BackgroundExpressionSampler[] geneBackgroundExpressionSamplers, SamplingRecords samplingRecords,
                               int threadNumber, int samplingTime, int burnIn, String tmpDir, String outputFile) {
        this.ipReads = ipReads;
        this.inputReads = inputReads;
        this.individualNumber = ipReads.length;
        this.peakNumber = ipReads[0].length;
        this.methylationLevelSampler = methLevelSampler;
        this.expansionEffectSampler = expansionEffectSampler;
        this.ipOverdispersionSampler = ipOverdispersionSampler;
        this.inputOverdispersionSampler = inputOverdispersionSampler;
        this.geneBackgroundExpressionSamplers = geneBackgroundExpressionSamplers;
        this.samplingRecords = samplingRecords;
        this.threadNumber = threadNumber;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.tmpDir = tmpDir;
        this.outputFile = outputFile;
    }

    /**
     * calculate IP and INPUT size factor for each individual
     */
    public void calcSampleSizeFactors() {
        // shape 1 × individualNumber
        double[] sampleIPSizeFactors = new double[this.individualNumber];
        double[] sampleINPUTSizeFactors = new double[this.individualNumber];
        SizeFactor sizeFactor;
        double ipFactor, inputFactor;
        for (int i=0; i<this.individualNumber; i++) {
            int[] individualIPReads = this.ipReads[i];
            int[] individualINPUTReads = this.inputReads[i];
            sizeFactor = new SizeFactor(individualIPReads, individualINPUTReads);
            ipFactor = sizeFactor.getSizeFactor(false);
            inputFactor = sizeFactor.getSizeFactor(true);
            sampleIPSizeFactors[i] = ipFactor;
            sampleINPUTSizeFactors[i] = inputFactor;
        }

        this.samplingRecords.setInputSizeFactors(sampleINPUTSizeFactors);
        this.samplingRecords.setIpSizeFactors(sampleIPSizeFactors);
    }

    /**
     * sampling process
     */
    public void iteration() {
        this.executorService = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
        this.renewRecordFiles();
        for (int t=1; t<this.samplingTime; t++) {
            this.inputOverdispersionSampling(this.samplingRecords, t);
            this.ipOverdispersionSampling(this.samplingRecords, t);
            this.expansionEffectSampling(this.samplingRecords, t);
            this.backGroundExpressionSampling(this.samplingRecords, t);
            this.methylationLevelSampling(this.samplingRecords, t);
            if (t<this.burnIn)
                this.renewRecordFiles();
        }
    }

    /**
     * return quantification result
     * @return quantification result
     */
    public HashMap<Integer, String> getQuantifyResult() {
        return this.quantifyResult;
    }

    /**
     * new sampling iteration of INPUT overdispersion for genes
     */
    private void inputOverdispersionSampling(SamplingRecords samplingRecords, int time) {

        double prevINPUTOverdispersion, curINPUTOverdispersion, curIPOverdispersion, curExpansion, prevPosteriorProba, curPosteriorProba;
        double[] curMethylationLevels, backgroundExpressions;  // shape 1 × peakNumber
        boolean samplingRes;

        prevINPUTOverdispersion = samplingRecords.getInputOverdispersionValue(time-1);
        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(prevINPUTOverdispersion);
        curIPOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);
        curExpansion = this.samplingRecords.getExpansionEffectsValue(time-1);
        curMethylationLevels = samplingRecords.getMethylationLevels();
        backgroundExpressions = samplingRecords.getBackgroundExpressions();

        prevPosteriorProba = this.calcLogLikelihood(curIPOverdispersion, prevINPUTOverdispersion, curExpansion,
                                                    curMethylationLevels, backgroundExpressions) +
                             this.calcLogPrior(curIPOverdispersion, prevINPUTOverdispersion, curExpansion,
                                               curMethylationLevels, backgroundExpressions);

        curPosteriorProba = this.calcLogLikelihood(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                                   curMethylationLevels, backgroundExpressions) +
                            this.calcLogPrior(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                              curMethylationLevels, backgroundExpressions);

        samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            samplingRecords.setInputOverdispersionValue(curINPUTOverdispersion, time);
        } else { // otherwise, rejected
            samplingRecords.setInputOverdispersionValue(prevINPUTOverdispersion, time);
        }
    }

    /**
     * new sampling iteration of INPUT overdispersion for genes
     */
    private void ipOverdispersionSampling(SamplingRecords samplingRecords, int time) {

        double curINPUTOverdispersion, prevIPOverdispersion, curIPOverdispersion, curExpansion, prevPosteriorProba, curPosteriorProba;
        double[] curMethylationLevels, backgroundExpressions;  // shape 1 × peakNumber
        boolean samplingRes;

        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        prevIPOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);
        curIPOverdispersion = this.inputOverdispersionSampler.randomSample(prevIPOverdispersion);
        curExpansion = this.samplingRecords.getExpansionEffectsValue(time-1);
        curMethylationLevels = samplingRecords.getMethylationLevels();
        backgroundExpressions = samplingRecords.getBackgroundExpressions();

        prevPosteriorProba = this.calcLogLikelihood(prevIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                                    curMethylationLevels, backgroundExpressions) +
                             this.calcLogPrior(prevIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                               curMethylationLevels, backgroundExpressions);

        curPosteriorProba = this.calcLogLikelihood(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                                   curMethylationLevels, backgroundExpressions) +
                            this.calcLogPrior(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                              curMethylationLevels, backgroundExpressions);

        samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            samplingRecords.setIpOverdispersionValue(curIPOverdispersion, time);
        } else { // otherwise, rejected
            samplingRecords.setIpOverdispersionValue(prevIPOverdispersion, time);
        }
    }

    /**
     * new sampling iteration of expansion effect of antigen
     */
    private void expansionEffectSampling(SamplingRecords samplingRecords, int time) {

        double curINPUTOverdispersion, curIPOverdispersion, prevExpansion, curExpansion, prevPosteriorProba, curPosteriorProba;
        double[] curMethylationLevels, backgroundExpressions;  // shape 1 × peakNumber
        boolean samplingRes;

        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        curIPOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        prevExpansion = this.samplingRecords.getExpansionEffectsValue(time-1);
        curExpansion = this.expansionEffectSampler.randomSample(prevExpansion);
        curMethylationLevels = samplingRecords.getMethylationLevels();
        backgroundExpressions = samplingRecords.getBackgroundExpressions();

        prevPosteriorProba = this.calcLogLikelihood(curIPOverdispersion, curINPUTOverdispersion, prevExpansion,
                                                    curMethylationLevels, backgroundExpressions) +
                             this.calcLogPrior(curIPOverdispersion, curINPUTOverdispersion, prevExpansion,
                                               curMethylationLevels, backgroundExpressions);

        curPosteriorProba = this.calcLogLikelihood(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                                   curMethylationLevels, backgroundExpressions) +
                            this.calcLogPrior(curIPOverdispersion, curINPUTOverdispersion, curExpansion,
                                              curMethylationLevels, backgroundExpressions);

        samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            samplingRecords.setExpansionEffectsValue(curExpansion, time);
        } else { // otherwise, rejected
            samplingRecords.setExpansionEffectsValue(prevExpansion, time);
        }
    }

    /**
     * new sampling iteration of background expression for each peak
     */
    private void backGroundExpressionSampling(SamplingRecords samplingRecords, int time) {

        double curINPUTOverdispersion, curIPOverdispersion, curExpansion;
        double[] curMethylationLevels, prevBackgroundExpressions, newSamplingResult;  // shape 1 × peakNumber

        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        curIPOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        curExpansion = this.samplingRecords.getExpansionEffectsValue(time);
        double[] sampleIPSizeFactors = this.samplingRecords.getIpSizeFactors();
        double[] sampleINPUTSizeFactors = this.samplingRecords.getInputSizeFactors();
        curMethylationLevels = samplingRecords.getMethylationLevels();
        prevBackgroundExpressions = samplingRecords.getBackgroundExpressions();
        newSamplingResult = new double[this.peakNumber];

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        SamplingTask task = ((prevBkgExp, curBkgExp, idx) -> () -> {
            double prevPosterior, curPosterior;
            boolean samplingRes;
            try {
                BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[idx];
                prevPosterior = this.calcPeakLogLikelihood(idx, curMethylationLevels[idx], prevBkgExp, curExpansion,
                                                           curIPOverdispersion, curINPUTOverdispersion,
                                                           sampleIPSizeFactors, sampleINPUTSizeFactors) +
                                sampler.getLogDensity(prevBkgExp);

                curPosterior = this.calcPeakLogLikelihood(idx, curMethylationLevels[idx], curBkgExp, curExpansion,
                                                          curIPOverdispersion, curINPUTOverdispersion,
                                                          sampleIPSizeFactors, sampleINPUTSizeFactors) +
                               sampler.getLogDensity(curBkgExp);

                samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosterior, prevPosterior, true);
                this.lock.lock();
                if (samplingRes)
                    newSamplingResult[idx] = curBkgExp;
                else // otherwise, rejected
                    newSamplingResult[idx] = prevBkgExp;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                this.lock.unlock();
                countDownLatch.countDown();
            }
        });

        double prevExpression, curExpression;
        BackgroundExpressionSampler sampler;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            prevExpression = prevBackgroundExpressions[peakIdx];
            sampler = this.geneBackgroundExpressionSamplers[peakIdx];

            curExpression = sampler.randomSample(prevExpression);
            runnable = task.getTask(prevExpression, curExpression, peakIdx);
            this.executorService.submit(runnable);
        }

        try {
            countDownLatch.await();
            samplingRecords.setBackgroundExpressions(newSamplingResult);
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdownNow();
        }
    }

    /**
     * new sampling iteration of methylation level for each peak
     */
    private void methylationLevelSampling(SamplingRecords samplingRecords, int time) {

        double curINPUTOverdispersion, curIPOverdispersion, curExpansion;
        double[] prevMethylationLevels, curBackgroundExpressions, newSamplingResult;  // shape 1 × peakNumber

        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        curIPOverdispersion = this.inputOverdispersionSampler.randomSample(time);
        curExpansion = this.samplingRecords.getExpansionEffectsValue(time);
        double[] sampleIPSizeFactors = this.samplingRecords.getIpSizeFactors();
        double[] sampleINPUTSizeFactors = this.samplingRecords.getInputSizeFactors();
        prevMethylationLevels = samplingRecords.getMethylationLevels();
        curBackgroundExpressions = samplingRecords.getBackgroundExpressions();
        newSamplingResult = new double[this.peakNumber];

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        SamplingTask task = ((prevMeth, curMeth, idx) -> () -> {
            double prevPosterior, curPosterior;
            boolean samplingRes;
            try {
                prevPosterior = this.calcPeakLogLikelihood(idx, prevMeth, curBackgroundExpressions[idx], curExpansion,
                                                           curIPOverdispersion, curINPUTOverdispersion,
                                                           sampleIPSizeFactors, sampleINPUTSizeFactors) +
                                this.methylationLevelSampler.getLogDensity(prevMeth);

                curPosterior = this.calcPeakLogLikelihood(idx, curMeth, curBackgroundExpressions[idx], curExpansion,
                                                          curIPOverdispersion, curINPUTOverdispersion,
                                                          sampleIPSizeFactors, sampleINPUTSizeFactors) +
                               this.methylationLevelSampler.getLogDensity(curMeth);

                samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosterior, prevPosterior, true);
                this.lock.lock();
                if (samplingRes)
                    newSamplingResult[idx] = curMeth;
                else // otherwise, rejected
                    newSamplingResult[idx] = prevMeth;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                this.lock.unlock();
                countDownLatch.countDown();
            }
        });

        double prevMethylationLevel, curMethylationLevel;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            prevMethylationLevel = prevMethylationLevels[peakIdx];
            curMethylationLevel = this.methylationLevelSampler.randomSample(prevMethylationLevel);

            runnable = task.getTask(prevMethylationLevel, curMethylationLevel, peakIdx);
            this.executorService.submit(runnable);
        }

        try {
            countDownLatch.await();
            samplingRecords.setMethylationLevels(newSamplingResult);
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdownNow();
        }
    }

    /**
     * log-likelihood
     * @param ipOverdispersion overdispersion of IP data
     * @param inputOverdispersion overdispersion of INPUT data
     * @param expansion expansion effect of antigen
     * @param methylationLevels methylation levels, shape 1 × peakNumber
     * @param backgroundExpressions background expressions, shape 1 × peakNumber
     * @return log likelihood
     */
    private double calcLogLikelihood(double ipOverdispersion, double inputOverdispersion, double expansion,
                                     double[] methylationLevels, double[] backgroundExpressions) {
        // shape 1 × individualNumber
        double[] sampleIPSizeFactor = this.samplingRecords.getIpSizeFactors();
        double[] sampleINPUTSizeFactor = this.samplingRecords.getInputSizeFactors();

        double expression, methylationLevel, proba = 0;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            expression = backgroundExpressions[peakIdx];
            methylationLevel = methylationLevels[peakIdx];
            proba += this.calcPeakLogLikelihood(peakIdx, methylationLevel, expression, expansion,
                                                ipOverdispersion, inputOverdispersion, sampleIPSizeFactor, sampleINPUTSizeFactor);
        }

        return proba;
    }

    /**
     * log-likelihood for a single peak region
     * @param peakIdx peak index
     * @param methylationLevel methylation level
     * @param expression expression
     * @param expansion expansion
     * @param ipOverdispersion IP data overdispersion
     * @param inputOverdispersion INPUT data overdispersion
     * @param sampleIPSizeFactor IP size factors for each sample
     * @param sampleINPUTSizeFactor INPUT size factors for each sample
     * @return log-likelihood for a single peak region
     */
    private double calcPeakLogLikelihood(int peakIdx, double methylationLevel, double expression, double expansion,
                                         double ipOverdispersion, double inputOverdispersion,
                                         double[] sampleIPSizeFactor, double[] sampleINPUTSizeFactor) {
        double ipFactor, inputFactor, proba = 0;
        // shape 1 × individualNumber
        int[] ipReadsCount = new int[this.individualNumber], inputReadsCount = new int[this.individualNumber];
        double[] ipExpectation = new double[this.individualNumber], inputExpectation = new double[this.individualNumber];

        for (int sampleIdx=0; sampleIdx<this.individualNumber; sampleIdx++) {
            ipReadsCount[sampleIdx] = this.ipReads[sampleIdx][peakIdx];
            inputReadsCount[sampleIdx] = this.inputReads[sampleIdx][peakIdx];
            ipFactor = sampleIPSizeFactor[sampleIdx];
            inputFactor = sampleINPUTSizeFactor[sampleIdx];

            ipExpectation[sampleIdx] = this.calcReadsExpectation(ipFactor, expression, methylationLevel, expansion);
            inputExpectation[sampleIdx] = this.calcReadsExpectation(inputFactor, expression, 1, 1);
        }
        proba += ProbabilityCalculator.logNegativeProbability(ipReadsCount, ipExpectation, ipOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(inputReadsCount, inputExpectation, inputOverdispersion);
        ipReadsCount = null;
        inputReadsCount = null;
        ipExpectation = null;
        inputExpectation = null;

        return proba;
    }

    /**
     * calculate IP and INPUT reads count expectation
     * @param sizeFactor size factor
     * @param expression background expression
     * @param methylationLevel methylation level
     * @param expansion expansion effect of antigen
     * @return reads count expectation
     */
    private double calcReadsExpectation(double sizeFactor, double expression, double methylationLevel, double expansion) {
        return sizeFactor * expansion * expression * methylationLevel;
    }

    /**
     * calculate log prior
     * @param ipOverdispersion overdispersion of IP data
     * @param inputOverdispersion overdispersion of INPUT data
     * @param expansion expansion effect of antigen
     * @param methylationLevels methylation levels, shape 1 × peakNumber
     * @param backgroundExpressions background expressions, shape 1 × peakNumber
     * @return log prior
     */
    private double calcLogPrior(double ipOverdispersion, double inputOverdispersion, double expansion,
                                double[] methylationLevels, double[] backgroundExpressions) {

        double proba = 0;
        double methylationLevel, expression;
        BackgroundExpressionSampler expressionSampler;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            methylationLevel = methylationLevels[peakIdx];
            expression = backgroundExpressions[peakIdx];
            expressionSampler = this.geneBackgroundExpressionSamplers[peakIdx];

            proba += this.methylationLevelSampler.getLogDensity(methylationLevel);
            proba += expressionSampler.getLogDensity(expression);
        }
        proba += this.ipOverdispersionSampler.getLogDensity(ipOverdispersion);
        proba += this.inputOverdispersionSampler.getLogDensity(inputOverdispersion);
        proba += this.expansionEffectSampler.getLogDensity(expansion);

        return proba;
    }

    /**
     * renew sampling files for each peak
     */
    private void renewRecordFiles() {
        double[] methylationLevels = this.samplingRecords.getMethylationLevels();
        double[] backgroundExpressions = this.samplingRecords.getBackgroundExpressions();

        double methylationValue, expressionValue;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            methylationValue = methylationLevels[peakIdx];
            expressionValue = backgroundExpressions[peakIdx];

            this.recordData(peakIdx, methylationValue, expressionValue);
        }
    }

    /**
     * write sampling result into file
     * @param peakIdx peak number
     * @param methylationValue methylation level
     * @param expressionValue expression value
     */
    private void recordData(int peakIdx, double methylationValue, double expressionValue) {
        File recordFile = new File(this.tmpDir, peakIdx + ".txt");
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new FileWriter(recordFile, true));
            bfw.write(String.join("\t", new String[] {String.valueOf(methylationValue),
                                                                String.valueOf(expressionValue)}));
            bfw.newLine();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * quantify methylation level and expression for each peak
     */
    public void quantify() {
        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        this.quantifyResult = new HashMap<>();
        QuantifyTask task = ((idx -> () -> {
            File targetFile = new File(this.tmpDir, idx+".txt");
            BufferedReader bfr = null;
            double[] methylationValues, expressionValues;
            boolean succ;
            try {
                bfr = new BufferedReader(new InputStreamReader(new FileInputStream(targetFile)));
                methylationValues = new double[this.samplingTime - this.burnIn];
                expressionValues = new double[this.samplingTime - this.burnIn];
                int lineNum = 0;
                String line = "";
                String[] info;
                double methylation, expression;
                while (line!=null) {
                    line = bfr.readLine();
                    if (line!=null) {
                        info = line.split("\t");
                        methylation = Double.valueOf(info[0]);
                        expression = Double.valueOf(info[1]);
                        methylationValues[lineNum] = methylation;
                        expressionValues[lineNum] = expression;
                    }
                }
                double methylationValue = this.calcMedian(methylationValues);
                double expressionValue = this.calcMedian(expressionValues);
                succ = targetFile.delete();
                if (!succ)
                    System.out.println("can not delete temporary file " + targetFile.getAbsolutePath());
                this.lock.lock();
                this.quantifyResult.put(idx, String.join("\t", new String[] {this.df.format(methylationValue),
                                                                                       this.df.format(expressionValue)}));
            } catch (IOException ie) {
                ie.printStackTrace();
            } finally {
                countDownLatch.countDown();
                this.lock.unlock();
                if (bfr!=null) {
                    try {
                        bfr.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                methylationValues = null;
                expressionValues = null;
            }
        }));
        // quantify methylation level and expression for each peak
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            runnable = task.getTask(peakIdx);
            this.executorService.submit(runnable);
        }
        try {
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
        } finally {
            this.executorService.shutdown();
        }

        double[] ipOverdispersions = this.samplingRecords.getIpOverdispersions();
        double[] inputOverdispersions = this.samplingRecords.getInputOverdispersions();
        double[] expansionEffects = this.samplingRecords.getExpansionEffects();
        this.ipOverdispersion = this.calcMedian(ipOverdispersions);
        this.inputOverdispersion = this.calcMedian(inputOverdispersions);
        this.expansionEffect = this.calcMedian(expansionEffects);
    }

    /**
     * calculate median values
     * @param values array
     * @return median
     */
    private double calcMedian(double[] values) {
        if (values.length == 1)
            return values[0];
        Arrays.sort(values);
        int medianIdx = values.length / 2;
        if (values.length % 2 == 0)
            return (values[medianIdx-1] + values[medianIdx]) * 0.5;
        else
            return values[medianIdx];
    }

    /**
     * output quantification result
     */
    public void output() {
        BufferedWriter bfw = null;
        String newLine, ipOverdispersion, inputOverdispersion, expansionEffect, record;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write("peakNumber\tmethylationLevel\texpression\tIPOverdispersion\tINPUTOverdispersion\texpansionEffect");
            bfw.newLine();
            ipOverdispersion = this.df.format(this.ipOverdispersion);
            inputOverdispersion = this.df.format(this.inputOverdispersion);
            expansionEffect = this.df.format(this.expansionEffect);

            for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
                record = this.quantifyResult.get(peakIdx);
                newLine = String.join("\t", new String[] {Integer.toString(peakIdx), record, ipOverdispersion, inputOverdispersion, expansionEffect});
                bfw.write(newLine);
                bfw.newLine();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw!=null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
