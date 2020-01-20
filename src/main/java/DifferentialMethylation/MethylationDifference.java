package DifferentialMethylation;

import LaplaceMetropolisEstimator.BestModelParams;
import LaplaceMetropolisEstimator.CovarianceMatrix;
import Quantification.BackgroundExpressionSampler;
import Quantification.ExpansionEffectSampler;
import Quantification.MethylationLevelSampler;
import Quantification.OverdispersionSampler;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.SizeFactor;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Judge whether the methylation level of two group of data, treatment and control, has difference.
 * H0: same methylation level
 * H1: difference in methylation level
 */
public class MethylationDifference {
    private int tretIndividualNumber, ctrlIndividualNumber, geneNumber, samplingTime, burnIn, threadNum;
    private int[][] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads;
    // distribution parameters
    private double ipOverdispersionShape, ipOverdispersionScale, inputOverdispersionShape, inputOverdispersionScale, expansionEffectParam1, expansionEffectParam2,
                   methylationLevelParam1, methylationLevelParam2, quantifiedExpansionEffect;
    // MH sampling components, shape geneNumber × samplingTime
    private double[] tretGeneBackgroundExpression, ctrlGeneBackgroundExpression;
    // size factors, shape 1 × individual number
    private double[] tretIPSizeFactors, tretINPUTSizeFactors, ctrlIPSizeFactors, ctrlINPUTSizeFactors;
    private String sameModelTmpDir, diffModelTmpDir, outputFile = null;
    private OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler;
    private ExpansionEffectSampler expansionEffectSampler;
    private MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler;
    private BackgroundExpressionSampler[] treatmentBackgroundExpressionSamplers, controlBackgroundExpressionSamplers;
    private HashMap<String, double[]> sameMethModelQuantify, diffMethModelQuantify;
    private ConcurrentHashMap<Integer, Double> bayesFactors;
    private ConcurrentHashMap<Integer, double[]> quantifyResult;
    private DecimalFormat df = new DecimalFormat("0.000000");

    /**
     * Constructor
     * @param treatmentIPReads treatment group IP reads count, shape individualNumber × geneNumber
     * @param treatmentINPUTReads treatment group INPUT reads count, shape individualNumber × geneNumber
     * @param controlIPReads control group IP reads count, shape individualNumber × geneNumber
     * @param controlINPUTReads control group INPUT reads count, shape individualNumber × geneNumber
     * @param ipOverdispersionShape IP data overdispersion Gamma distribution shape parameter
     * @param ipOverdispersionScale IP data overdispersion Gamma distribution scale parameter
     * @param inputOverdispersionShape INPUT data overdispersion Gamma distribution shape parameter
     * @param inputOverdispersionScale INPUT data overdispersion Gamma distribution scale parameter
     * @param methylationLevelParam1 methylation level Beta distribution shape parameter
     * @param methylationLevelParam2 methylation level Beta distribution shape parameter
     * @param samplingTime MH sampling time
     * @param burnIn MH burn-in time
     */
    public MethylationDifference(int[][] treatmentIPReads, int[][] treatmentINPUTReads,
                                 int[][] controlIPReads, int[][] controlINPUTReads,
                                 double ipOverdispersionShape, double ipOverdispersionScale,
                                 double inputOverdispersionShape, double inputOverdispersionScale,
                                 double expansionEffectParam1, double expansionEffectParam2,
                                 double methylationLevelParam1, double methylationLevelParam2,
                                 int samplingTime, int burnIn, String outputFile, int threadNumber) {
        // checking input parameters, same replications in IP and INPUT sample
        assert treatmentINPUTReads.length == treatmentIPReads.length && controlINPUTReads.length == controlIPReads.length;
        // treatment and control group contains same gene number in IP and INPUT
        assert treatmentINPUTReads[0].length == treatmentIPReads[0].length && controlINPUTReads[0].length == controlIPReads[0].length;
        assert treatmentINPUTReads[0].length == controlINPUTReads[0].length;
        assert samplingTime > burnIn;

        this.controlINPUTReads = controlINPUTReads;
        this.controlIPReads = controlIPReads;
        this.treatmentINPUTReads = treatmentINPUTReads;
        this.treatmentIPReads = treatmentIPReads;
        this.tretIndividualNumber = treatmentINPUTReads.length;
        this.ctrlIndividualNumber = controlINPUTReads.length;
        this.geneNumber = treatmentINPUTReads[0].length;
        this.ipOverdispersionShape = ipOverdispersionShape;
        this.ipOverdispersionScale = ipOverdispersionScale;
        this.inputOverdispersionShape = inputOverdispersionShape;
        this.inputOverdispersionScale = inputOverdispersionScale;
        this.expansionEffectParam1 = expansionEffectParam1;
        this.expansionEffectParam2 = expansionEffectParam2;
        this.methylationLevelParam1 = methylationLevelParam1;
        this.methylationLevelParam2 = methylationLevelParam2;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNum = threadNumber;
        this.outputFile = outputFile;
        this.sameModelTmpDir = new File(new File(this.outputFile).getParent(), "sameModel").getAbsolutePath();
        this.diffModelTmpDir = new File(new File(this.outputFile).getParent(), "diffModel").getAbsolutePath();
    }

    /**
     * API function
     */
    public void runExecute() {
        this.initialize();
        this.selectModel();
        if (this.outputFile != null)
            this.output(this.outputFile);
    }

    /**
     * initialize for MH sampling
     */
    private void initialize() {
        // initialize priority distribution
        this.tretMethylationLevelSampler = new MethylationLevelSampler(this.methylationLevelParam1, this.methylationLevelParam2);
        this.ctrlMethylationLevelSampler = new MethylationLevelSampler(this.methylationLevelParam1, this.methylationLevelParam2);
        this.tretIPOverdispersionSampler = new OverdispersionSampler(this.ipOverdispersionShape, this.ipOverdispersionScale);
        this.tretINPUTOverdispersionSampler = new OverdispersionSampler(this.inputOverdispersionShape, this.inputOverdispersionScale);
        this.ctrlIPOverdispersionSampler = new OverdispersionSampler(this.ipOverdispersionShape, this.ipOverdispersionScale);
        this.ctrlINPUTOverdispersionSampler = new OverdispersionSampler(this.inputOverdispersionShape, this.inputOverdispersionScale);
        this.expansionEffectSampler = new ExpansionEffectSampler(this.expansionEffectParam1, this.expansionEffectParam2);

        // background expression of treatment group, size factor of treatment group IP and INPUT data
        BackgroundExpression be = new BackgroundExpression(this.treatmentIPReads, this.treatmentINPUTReads);
        double[] backgroundExpressionMean = be.geneBackgroundExp(); // background expression expectation, shape 1 × geneNumber
        double[] backgroundExpressionStd = be.geneExpressionStd(); // background expression standard deviation, shape 1 × geneNumber
        double scale, shape;
        this.tretGeneBackgroundExpression = new double[this.geneNumber];
        this.treatmentBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            // shape and scale parameters of log-normal distribution
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = backgroundExpressionStd[geneIdx];
            this.tretGeneBackgroundExpression[geneIdx] = backgroundExpressionMean[geneIdx];
            this.treatmentBackgroundExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }
        this.tretIPSizeFactors = be.getGlobalIPSizeFactor();
        this.tretINPUTSizeFactors = be.getGlobalINPUTSizeFactor();

        // background expression of control group, size factor of control group IP and INPUT data
        be = new BackgroundExpression(this.controlIPReads, this.controlINPUTReads);
        backgroundExpressionMean = be.geneBackgroundExp();
        backgroundExpressionStd = be.geneExpressionStd();
        this.ctrlGeneBackgroundExpression = new double[this.geneNumber];
        this.controlBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = backgroundExpressionStd[geneIdx];
            this.ctrlGeneBackgroundExpression[geneIdx] = backgroundExpressionMean[geneIdx];
            this.controlBackgroundExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }
        this.ctrlIPSizeFactors = be.getGlobalIPSizeFactor();
        this.ctrlINPUTSizeFactors = be.getGlobalINPUTSizeFactor();

        // initialize samples' size factors
        this.calcSampleSizeFactors();
    }

    /**
     * pseudo prior distribution for IP and INPUT Overdispersion, methylation level, and background expression for each gene
     */
    private void selectModel() {
        this.bayesFactors = new ConcurrentHashMap<>();
        this.quantifyResult = new ConcurrentHashMap<>();

        SameMethylationLevelModel sameMethylationLevelModel;
        DiffMethylationLevelModel quantifyExpansionEffectModel, diffMethylationLevelModel;
        ExecutorService executorService = Executors.newFixedThreadPool(this.threadNum);
        try {
            quantifyExpansionEffectModel = (DiffMethylationLevelModel) this.initialModel(false);
            // quantify expansion effect of treatment and control IP data
            this.quantifyExpansionEffect(quantifyExpansionEffectModel);

            sameMethylationLevelModel = (SameMethylationLevelModel) this.initialModel(true);
            sameMethylationLevelModel.setQuantifiedExpansionEffect(this.quantifiedExpansionEffect);
            diffMethylationLevelModel = (DiffMethylationLevelModel) this.initialModel(false);
            diffMethylationLevelModel.setQuantifiedExpansionEffect(this.quantifiedExpansionEffect);
            // quantify treatment and control group methylation level with specified expansion effect
            this.differentiateSampling(sameMethylationLevelModel, diffMethylationLevelModel);
            this.sameMethModelQuantify = this.quantify(sameMethylationLevelModel);
            this.diffMethModelQuantify = this.quantify(diffMethylationLevelModel);

            CountDownLatch countDownLatch = new CountDownLatch(this.geneNumber);

            // calculate Bayes factors for each peak
            BayesFactorTask taskMaker = (idx -> () -> {
                try {
                    double bayesFactor = this.calcBayesFactorLaplaceMetropolis(sameMethylationLevelModel, diffMethylationLevelModel, idx);
                    double tretMeth, ctrlMeth, tretExp, ctrlExp, tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion;
                    ;
                    if (bayesFactor - 3 > 0.00001) {
                        tretMeth = this.diffMethModelQuantify.get("tretMeth")[idx];
                        ctrlMeth = this.diffMethModelQuantify.get("ctrlMeth")[idx];
                        tretExp = this.diffMethModelQuantify.get("tretExp")[idx];
                        ctrlExp = this.diffMethModelQuantify.get("ctrlExp")[idx];
                        tretIPOverdispersion = this.diffMethModelQuantify.get("overdispersions")[0];
                        ctrlIPOverdispersion = this.diffMethModelQuantify.get("overdispersions")[1];
                        tretINPUTOverdispersion = this.diffMethModelQuantify.get("overdispersions")[2];
                        ctrlINPUTOverdispersion = this.diffMethModelQuantify.get("overdispersions")[3];
                    } else {
                        tretMeth = this.sameMethModelQuantify.get("tretMeth")[idx];
                        ctrlMeth = this.sameMethModelQuantify.get("ctrlMeth")[idx];
                        tretExp = this.sameMethModelQuantify.get("tretExp")[idx];
                        ctrlExp = this.sameMethModelQuantify.get("ctrlExp")[idx];
                        tretIPOverdispersion = this.sameMethModelQuantify.get("overdispersions")[0];
                        ctrlIPOverdispersion = this.sameMethModelQuantify.get("overdispersions")[1];
                        tretINPUTOverdispersion = this.sameMethModelQuantify.get("overdispersions")[2];
                        ctrlINPUTOverdispersion = this.sameMethModelQuantify.get("overdispersions")[3];
                    }
                    this.bayesFactors.put(idx, bayesFactor);
                    this.quantifyResult.put(idx, new double[] {tretMeth, ctrlMeth, tretExp, ctrlExp, tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion});
                } catch (Exception e) {
                    e.printStackTrace();
                } finally {
                    countDownLatch.countDown();
                }
            });

            Runnable runnable;
            for (int peakIdx=0; peakIdx<this.geneNumber; peakIdx++) {
                runnable = taskMaker.createTask(peakIdx);
                executorService.submit(runnable);
            }
            countDownLatch.await();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            executorService.shutdown();
        }
    }

    /**
     * initialize same methylation or diff methylation model
     * @param sameMeth true, if same model; otherwise false
     * @return model
     */
    private ModelSelection initialModel(boolean sameMeth) {
        ModelSelection model;
        if (sameMeth)
            model = new SameMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                  this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                  this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                  this.treatmentBackgroundExpressionSamplers, this.controlBackgroundExpressionSamplers,
                                                  this.expansionEffectSampler, this.samplingTime, this.burnIn, this.threadNum, this.sameModelTmpDir);
        else
            model = new DiffMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                  this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                  this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                  this.treatmentBackgroundExpressionSamplers, this.controlBackgroundExpressionSamplers,
                                                  this.expansionEffectSampler, this.samplingTime, this.burnIn, this.threadNum, this.diffModelTmpDir);
        // initial model samples' size factors
        model.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors);
        // initial model samples' observed reads counts
        model.setTreatmentControlGeneReads(this.treatmentIPReads, this.treatmentINPUTReads, this.controlIPReads, this.controlINPUTReads);
        // set sampling list
        double[] tretIPOverdispersion = new double[this.samplingTime], ctrlIPOverdispersion = new double[this.samplingTime],
                 tretINPUTOverdispersion = new double[this.samplingTime], ctrlINPUTOverdispersion = new double[this.samplingTime],
                 expansionEffect = new double[this.samplingTime], tretMethLevel = new double[this.geneNumber], ctrlMethLevel = new double[this.geneNumber];
        tretIPOverdispersion[0] = 0.1;
        ctrlIPOverdispersion[0] = 0.1;
        tretINPUTOverdispersion[0] = 0.1;
        ctrlINPUTOverdispersion[0] = 0.1;
        expansionEffect[0] = 2;
        tretMethLevel[0] = 0.1;
        ctrlMethLevel[0] = 0.1;

        double[] tretBkgExp = new double[this.geneNumber];
        System.arraycopy(this.tretGeneBackgroundExpression, 0, tretBkgExp, 0, this.geneNumber);
        double[] ctrlBkgExp = new double[this.geneNumber];
        System.arraycopy(this.ctrlGeneBackgroundExpression, 0, ctrlBkgExp, 0, this.geneNumber);
        model.setSamplingList(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                              expansionEffect, tretBkgExp, ctrlBkgExp, tretMethLevel, ctrlMethLevel);

        return model;
    }

    /**
     * quantify expansion effect of treatment and control data
     * @param diffMethylationLevelModel different methylation level model
     */
    private void quantifyExpansionEffect(DiffMethylationLevelModel diffMethylationLevelModel) {
        double diffMethFullConditionProba;
        diffMethFullConditionProba = diffMethylationLevelModel.paramsLogFullConditionalProbability(0);
        for (int t=1; t<this.samplingTime; t++) {
            diffMethFullConditionProba = diffMethylationLevelModel.iterate(t, diffMethFullConditionProba);
        }
        double[] expansionEffect = diffMethylationLevelModel.parameters.getExpansionEffect();
        this.quantifiedExpansionEffect = diffMethylationLevelModel.calcMedian(expansionEffect);
    }

    /**
     * sampling process
     */
    private void differentiateSampling(SameMethylationLevelModel sameMethylationLevelModel, DiffMethylationLevelModel diffMethylationLevelModel) {
        double sameMethFullConditionProba, diffMethFullConditionProba;
        // iterate with specific expansion effect, record model initial value
        sameMethylationLevelModel.renewRecordFiles();
        diffMethylationLevelModel.renewRecordFiles();
        // calculate model initial posterior probability, log-likelihood + log-prior
        sameMethFullConditionProba = sameMethylationLevelModel.paramsLogFullConditionalProbability(0);
        diffMethFullConditionProba = diffMethylationLevelModel.paramsLogFullConditionalProbability(0);
        for (int t=1; t<this.samplingTime; t++) {
            sameMethFullConditionProba = sameMethylationLevelModel.iterate(t, sameMethFullConditionProba);
            diffMethFullConditionProba = diffMethylationLevelModel.iterate(t, diffMethFullConditionProba);
            // record sampling result of this iteration
            sameMethylationLevelModel.renewRecordFiles();
            diffMethylationLevelModel.renewRecordFiles();
        }
    }

    /**
     * treatment and control group methylation, IP overdispersion, INPUT overdispersion, background expression quantification
     * use quantification model Quantification.MethylationQuantification
     * @param model model need to be quantified
     * @return double[] {treatment methylation, control methylation, treatment IP overdispersion, treatment INPUT overdispersion,
     *                   control IP overdispersion, control INPUT overdispersion, treatment background expression, control background expression}
     */
    private HashMap<String, double[]> quantify(ModelSelection model) {
        model.quantify();
        double[] overdispersions = model.parameters.getOverdispersions();
        // shape 1 × peakNumber
        double[] tretExpressions = model.parameters.getQuantifiedTretExpression();
        double[] ctrlExpressions = model.parameters.getQuantifiedCtrlExpression();
        double[] tretMethylations = model.parameters.getQuantifiedTretMethLevel();
        double[] ctrlMethylations = model.parameters.getQuantifiedCtrlMethLevel();

        HashMap<String, double[]> quantifyResult = new HashMap<>();
        quantifyResult.put("overdispersions", overdispersions);
        quantifyResult.put("tretExp", tretExpressions);
        quantifyResult.put("ctrlExp", ctrlExpressions);
        quantifyResult.put("tretMeth", tretMethylations);
        quantifyResult.put("ctrlMeth", ctrlMethylations);

        return quantifyResult;
    }

    /**
     * Deprecated
     * probability of variant methylation level between treatment and control
     * @return probability
     */
    @Deprecated
    private double differentiation(int[] modelSelectionResult) {
        int remainTime = this.samplingTime - this.burnIn;
        return  (double) (Arrays.stream(modelSelectionResult).skip(this.burnIn).sum()) / remainTime;
    }

    /**
     * estimate Bayes factor with Laplace method, need second-order derivative calculation
     * Laplace method estimate the marginal probability p(D|M_k),
     *              p(D|M1)
     *      BF = -------------  ->  log BF = p(D|M1) - p(D|M0)
     *              p(D|M0)
     *
     *  2*log(BF) is used as a evidential measure to compare the support provided by the observe data D for Model0 relative to Model1.
     */
    @Deprecated
    private double calcLogBayesFactorLaplaceApproximation(double[] sameMethModelQuantify, double[] diffMethModelQuantify,
                                                          ModelSelection sameMethModel, ModelSelection diffMethModel, int geneIdx) {
//        int[] tretIPReads = sameMethModel.treatmentIPReads;
//        int[] tretINPUTReads = sameMethModel.treatmentINPUTReads;
//        int[] ctrlIPReads = sameMethModel.controlIPReads;
//        int[] ctrlINPUTReads = sameMethModel.controlINPUTReads;
//        int[] tretIPNonPeak = sameMethModel.treatmentIPNonPeakReads;
//        int[] tretINPUTNonPeak = sameMethModel.treatmentINPUTNonPeakReads;
//        int[] ctrlIPNonPeak = sameMethModel.controlIPNonPeakReads;
//        int[] ctrlINPUTNonPeak = sameMethModel.controlINPUTNonPeakReads;
//        double tretPeakBkgExp = this.tretGeneBackgroundExpression[geneIdx];
//        double tretNonPeakBkgExp= this.tretNonPeakExpressionMean[geneIdx];
//        double ctrlPeakBkgExp = this.ctrlGeneBackgroundExpression[geneIdx];
//        double ctrlNonPeakBkgExp= this.ctrlNonPeakExpressionMean[geneIdx];
//        double model0LogLikelihood = this.modelLogLikelihood(sameMethModelQuantify, sameMethModel);
//        double model1LofLikelihood = this.modelLogLikelihood(diffMethModelQuantify, diffMethModel);
//        double model0LogPrior = this.modelLogPrior(sameMethModelQuantify, sameMethModel);
//        double model1LogPrior = this.modelLogPrior(diffMethModelQuantify, diffMethModel);
//        double[] tretSizeFactor = this.sampleSizeFactors(this.treatmentIPReads, this.treatmentINPUTReads);
//        double[] ctrlSizeFactor = this.sampleSizeFactors(this.controlIPReads, this.controlINPUTReads);
//        BayesFactorCalculation bfc = new BayesFactorCalculation(model0LogLikelihood, model1LofLikelihood, model0LogPrior, model1LogPrior,
//                                                                sameMethModelQuantify, diffMethModelQuantify, tretPeakBkgExp, ctrlPeakBkgExp, tretNonPeakBkgExp, ctrlNonPeakBkgExp,
//                                                                this.ipOverdispersionShape, this.methylationLevelParam1, this.methylationLevelParam2);
//        bfc.setReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
//                     tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak);
//        bfc.setSizeFactors(tretSizeFactor, ctrlSizeFactor);
//        return bfc.calcBayesFactor();
        return 0;
    }

    /**
     * calculate Bayes factor with Laplace-Metropolis estimator
     * @param sameMethModel same methylation level model
     * @param diffMethModel diff methylation level model
     * @param peakIdx peakIdx
     * @return Bayes factor
     */
    private double calcBayesFactorLaplaceMetropolis(ModelSelection sameMethModel, ModelSelection diffMethModel, int peakIdx) throws IOException {
        // shape 1 × samplingTime, get from model
        double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion;
        // shape 1 × samplingTime, get from record file
        double[][] recordData;
        double[] treatmentMethylationLevel, controlMethylationLevel, treatmentExpressions, controlExpressions;
        File sameModelTargetFile = new File(this.sameModelTmpDir, peakIdx+".txt"),
             diffModelTargetFile = new File(this.diffModelTmpDir, peakIdx+".txt");

        // same methylation level model
        treatmentIPOverdispersion = sameMethModel.parameters.getTretIPOverdispersion();
        treatmentINPUTOverdispersion = sameMethModel.parameters.getTretINPUTOverdispersion();
        controlIPOverdispersion = sameMethModel.parameters.getCtrlIPOverdispersion();
        controlINPUTOverdispersion = sameMethModel.parameters.getCtrlINPUTOverdispersion();
        recordData = this.loadRecordData(sameModelTargetFile);
        treatmentMethylationLevel = recordData[0];
        controlMethylationLevel = recordData[1];
        treatmentExpressions = recordData[2];
        controlExpressions = recordData[3];

        double model0Determinate = this.calcCovarianceMatrixDeterminate(treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                                                                        treatmentMethylationLevel, null, treatmentExpressions, controlExpressions);

        double model0MaximumPosterior = this.calcPeakMaximumPosterior(treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                                                                      treatmentMethylationLevel, null, treatmentExpressions, controlExpressions, sameMethModel, peakIdx);

        // diff methylation level model
        treatmentIPOverdispersion = diffMethModel.parameters.getTretIPOverdispersion();
        treatmentINPUTOverdispersion = diffMethModel.parameters.getTretINPUTOverdispersion();
        controlIPOverdispersion = diffMethModel.parameters.getCtrlIPOverdispersion();
        controlINPUTOverdispersion = diffMethModel.parameters.getCtrlINPUTOverdispersion();
        recordData = this.loadRecordData(diffModelTargetFile);
        treatmentMethylationLevel = recordData[0];
        controlMethylationLevel = recordData[1];
        treatmentExpressions = recordData[2];
        controlExpressions = recordData[3];

        double model1Determinate = this.calcCovarianceMatrixDeterminate(treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                                                                        treatmentMethylationLevel, controlMethylationLevel, treatmentExpressions, controlExpressions);

        double model1MaximumPosterior = this.calcPeakMaximumPosterior(treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                                                                      treatmentMethylationLevel, controlMethylationLevel, treatmentExpressions, controlExpressions, diffMethModel, peakIdx);

        int model0ParamNum = 7, model1ParamNum = 8;

        // p(D|M=0) and p(D|M=1)
        double model0Proba = Math.pow(2*Math.PI, model0ParamNum*0.5) * Math.pow(model0Determinate, 0.5) * model0MaximumPosterior;
        double model1Proba = Math.pow(2*Math.PI, model1ParamNum*0.5) * Math.pow(model1Determinate, 0.5) * model1MaximumPosterior;

        return model1Proba / model0Proba;
    }

    /**
     * calculate model covariance matrix determinate
     * @param treatmentIPOverdispersion treatment IP Overdispersion, shape 1 × samplingTime
     * @param treatmentINPUTOverdispersion treatment INPUT Overdispersion, shape 1 × samplingTime
     * @param controlIPOverdispersion control IP Overdispersion, shape 1 × samplingTime
     * @param controlINPUTOverdispersion control INPUT Overdispersion, shape 1 × samplingTime
     * @param treatmentMethylationLevel treatment methylation level, shape 1 × samplingTime
     * @param controlMethylationLevel control methylation level, shape 1 × samplingTime
     * @param treatmentExpressions treatment expression, shape 1 × samplingTime
     * @param controlExpressions control expression, shape 1 × samplingTime
     * @return covariance matrix determinate
     */
    private double calcCovarianceMatrixDeterminate(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                                                   double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                                                   double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                                                   double[] treatmentExpressions, double[] controlExpressions) {
        CovarianceMatrix modelCovariance = new CovarianceMatrix(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                                 controlIPOverdispersion, controlINPUTOverdispersion,
                                                                 treatmentMethylationLevel, controlMethylationLevel,
                                                                 treatmentExpressions, controlExpressions, this.burnIn);
        double[][] modelCovarianceMatrix = modelCovariance.getCovarianceMatrix();
        RealMatrix modelCovMatrix = new Array2DRowRealMatrix(modelCovarianceMatrix);

        return new LUDecomposition(modelCovMatrix).getDeterminant();
    }

    /**
     * calculate maximum posterior probability for a single peak
     * @param treatmentIPOverdispersion treatment IP Overdispersion, shape 1 × samplingTime
     * @param treatmentINPUTOverdispersion treatment INPUT Overdispersion, shape 1 × samplingTime
     * @param controlIPOverdispersion control IP Overdispersion, shape 1 × samplingTime
     * @param controlINPUTOverdispersion control INPUT Overdispersion, shape 1 × samplingTime
     * @param treatmentMethylationLevel treatment methylation level, shape 1 × samplingTime
     * @param controlMethylationLevel control methylation level, shape 1 × samplingTime
     * @param treatmentExpressions treatment expression, shape 1 × samplingTime
     * @param controlExpressions control expression, shape 1 × samplingTime
     * @param model model
     * @param peakIdx peak index
     * @return maximum posterior
     */
    private double calcPeakMaximumPosterior(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                                            double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                                            double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                                            double[] treatmentExpressions, double[] controlExpressions,
                                            ModelSelection model, int peakIdx) {
        BestModelParams modelBestParams = new BestModelParams(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                              controlIPOverdispersion, controlINPUTOverdispersion,
                                                              treatmentMethylationLevel, controlMethylationLevel,
                                                              treatmentExpressions, controlExpressions,
                                                              model, this.quantifiedExpansionEffect, this.burnIn, peakIdx);
        modelBestParams.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors,
                                       this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors);

        return modelBestParams.getMaximumPosterior();
    }

    /**
     * load sampling result of a peak
     * @param targetFile record file
     * @return record data
     */
    private double[][] loadRecordData(File targetFile) throws IOException {
        BufferedReader bfr = null;
        double[] tretMethList = new double[this.samplingTime], ctrlMethList = new double[this.samplingTime],
                 tretExpList = new double[this.samplingTime], ctrlExpList = new double[this.samplingTime];
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(targetFile)));
            String line = "";
            int lineNum = 0;
            double tretMeth, ctrlMeth, tretExp, ctrlExp;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    info = line.split("\t");
                    tretMeth = Double.valueOf(info[0]);
                    ctrlMeth = Double.valueOf(info[1]);
                    tretExp = Double.valueOf(info[2]);
                    ctrlExp = Double.valueOf(info[3]);

                    tretMethList[lineNum] = tretMeth;
                    tretExpList[lineNum] = tretExp;
                    ctrlMethList[lineNum] = ctrlMeth;
                    ctrlExpList[lineNum] = ctrlExp;
                    lineNum++;
                }
            }
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return new double[][] {tretMethList, ctrlMethList, tretExpList, ctrlExpList};
    }

    /**
     * size factors of each individual
     */
    private void calcSampleSizeFactors() {
        // shape 1 × individualNumber
        this.tretIPSizeFactors = new double[this.tretIndividualNumber];
        this.tretINPUTSizeFactors = new double[this.tretIndividualNumber];
        this.ctrlIPSizeFactors = new double[this.ctrlIndividualNumber];
        this.ctrlINPUTSizeFactors = new double[this.ctrlIndividualNumber];
        SizeFactor sizeFactor;
        double ipFactor, inputFactor;
        int[] individualIPReads, individualINPUTReads;
        for (int i=0; i<this.tretIndividualNumber; i++) {
            individualIPReads = this.treatmentIPReads[i];
            individualINPUTReads = this.treatmentINPUTReads[i];
            sizeFactor = new SizeFactor(individualIPReads, individualINPUTReads);
            ipFactor = sizeFactor.getSizeFactor(false);
            inputFactor = sizeFactor.getSizeFactor(true);
            this.tretIPSizeFactors[i] = ipFactor;
            this.tretINPUTSizeFactors[i] = inputFactor;
        }

        for (int i=0; i<this.ctrlIndividualNumber; i++) {
            individualIPReads = this.controlIPReads[i];
            individualINPUTReads = this.controlINPUTReads[i];
            sizeFactor = new SizeFactor(individualIPReads, individualINPUTReads);
            ipFactor = sizeFactor.getSizeFactor(false);
            inputFactor = sizeFactor.getSizeFactor(true);
            this.ctrlIPSizeFactors[i] = ipFactor;
            this.ctrlINPUTSizeFactors[i] = inputFactor;
        }
        // only for simulation data test
        // return new double[] {1, 1, 1};
    }

    /**
     * output test result into file
     * @param fileName output filename
     */
    private void output(String fileName) {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(fileName))));
            ArrayList<Integer> sortedList = new ArrayList<>(this.bayesFactors.keySet().stream().sorted(((o1, o2) -> o1-o2)).collect(Collectors.toList()));
            bfw.write("#geneIdx\tBF\ttretMeth\tctrlMeth\ttretBkgExp\tctrlBkgExp\ttretIPOverdispersion\ttretINPUTOverdispersion\tctrlIPOverdispersion\tctrlINPUTOverdispersion");//\ttretNonPeakExp\tctrlNonPeakExp
            bfw.newLine();
            String line, records, bf;
            double[] quantify;
            for (Integer geneIdx: sortedList) {
                quantify = this.quantifyResult.get(geneIdx);
                records = Arrays.stream(quantify).mapToObj(x -> (this.df.format((Double) x))).collect(Collectors.joining("\t"));
                bf = this.df.format(this.bayesFactors.get(geneIdx));
                line = String.join("\t",new String[] {Integer.toString(geneIdx+1), bf, records});
                bfw.write(line);
                bfw.newLine();
            }
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
}
