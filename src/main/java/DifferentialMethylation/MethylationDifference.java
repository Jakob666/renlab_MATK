package DifferentialMethylation;

import LaplaceMetropolisEstimator.BestModelParams;
import LaplaceMetropolisEstimator.CovarianceMatrix;
import Quantification.BackgroundExpressionSampler;
import Quantification.MethylationLevelSampler;
import Quantification.NonSpecificEnrichmentSampler;
import Quantification.OverdispersionSampler;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.ReadsExpectation;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Judge whether the methylation level of two group of data, treatment and control, has difference.
 * H0: same methylation level
 * H1: difference in methylation level
 */
public class MethylationDifference {
    private int tretIndividualNumber, ctrlIndividualNumber, geneNumber, samplingTime, burnIn, threadNum;
    private int[][] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads,
                    treatmentIPNonPeak, treatmentINPUTNonPeak, controlIPNonPeak, controlINPUTNonPeak;
    // distribution parameters
    private double ipOverdispersionShape, ipOverdispersionScale, inputOverdispersionShape, inputOverdispersionScale,
                   nonspecificEnrichParam1, nonspecificEnrichParam2, methylationLevelParam1, methylationLevelParam2, sameMethLevelModelPriorProba, diffMethLevelModelPriorProba;
    private double[][] tretIPReadsExpectation, tretINPUTReadsExpectation, ctrlIPReadsExpectation, ctrlINPUTReadsExpectation,
                       tretIPNonPeakExpectation, tretINPUTNonPeakExpectation, ctrlIPNonPeakExpectation, ctrlINPUTNonPeakExpectation;
    // MH sampling components, shape geneNumber × samplingTime
    private double[] tretGeneBackgroundExpression, ctrlGeneBackgroundExpression, tretNonPeakExpressionMean, ctrlNonPeakExpressionMean;
    // size factors, shape 1 × individual number
    private double[] tretIPSizeFactors, tretINPUTSizeFactors, tretIPNonPeakSizeFactors, tretINPUTNonPeakSizeFactors,
                     ctrlIPSizeFactors, ctrlINPUTSizeFactors, ctrlIPNonPeakSizeFactors, ctrlINPUTNonPeakSizeFactors;
    private String outputFile = null;
    private OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler;
    private NonSpecificEnrichmentSampler nonspecificEnrichSampler;
    private MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler;
    private BackgroundExpressionSampler[] treatmentBackgroundExpressionSamplers, controlBackgroundExpressionSamplers,
                                          treatmentNonPeakExpressionSamplers, controlNonPeakExpressionSamplers;
    private ModelSelector modelSelector;
    private ConcurrentHashMap<Integer, Double> bayesFactors;
    private ConcurrentHashMap<Integer, double[]> quantifyResult;
    private DecimalFormat df = new DecimalFormat("0.0000");

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
    public MethylationDifference(int[][] treatmentIPReads, int[][] treatmentINPUTReads, int[][] treatmentIPNonPeak, int[][] treatmentINPUTNonPeak,
                                 int[][] controlIPReads, int[][] controlINPUTReads, int[][] controlIPNonPeak, int[][] controlINPUTNonPeak,
                                 double ipOverdispersionShape, double ipOverdispersionScale,
                                 double inputOverdispersionShape, double inputOverdispersionScale,
                                 double nonspecificEnrichParam1, double nonspecificEnrichParam2,
                                 double methylationLevelParam1, double methylationLevelParam2,
                                 double sameMethLevelModelPriorProba, double diffMethLevelModelPriorProba,
                                 int samplingTime, int burnIn, String outputFile, int threadNumber) {
        // checking input parameters, same replications in IP and INPUT sample
        assert Math.abs(sameMethLevelModelPriorProba+diffMethLevelModelPriorProba-1) < 0.00001;
        assert treatmentINPUTReads.length == treatmentIPReads.length && controlINPUTReads.length == controlIPReads.length;
        // treatment and control group contains same gene number in IP and INPUT
        assert treatmentINPUTReads[0].length == treatmentIPReads[0].length && controlINPUTReads[0].length == controlIPReads[0].length;
        assert treatmentINPUTReads[0].length == controlINPUTReads[0].length;
        assert samplingTime > burnIn;

        this.controlINPUTReads = controlINPUTReads;
        this.controlIPReads = controlIPReads;
        this.treatmentINPUTReads = treatmentINPUTReads;
        this.treatmentIPReads = treatmentIPReads;
        this.controlIPNonPeak = controlIPNonPeak;
        this.controlINPUTNonPeak = controlINPUTNonPeak;
        this.treatmentIPNonPeak = treatmentIPNonPeak;
        this.treatmentINPUTNonPeak = treatmentINPUTNonPeak;
        this.tretIndividualNumber = treatmentINPUTReads.length;
        this.ctrlIndividualNumber = controlINPUTReads.length;
        this.geneNumber = treatmentINPUTReads[0].length;
        this.ipOverdispersionShape = ipOverdispersionShape;
        this.ipOverdispersionScale = ipOverdispersionScale;
        this.inputOverdispersionShape = inputOverdispersionShape;
        this.inputOverdispersionScale = inputOverdispersionScale;
        this.nonspecificEnrichParam1 = nonspecificEnrichParam1;
        this.nonspecificEnrichParam2 = nonspecificEnrichParam2;
        this.methylationLevelParam1 = methylationLevelParam1;
        this.methylationLevelParam2 = methylationLevelParam2;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNum = threadNumber;
        this.modelSelector = new ModelSelector();
        this.sameMethLevelModelPriorProba = sameMethLevelModelPriorProba;
        this.diffMethLevelModelPriorProba = diffMethLevelModelPriorProba;
        this.outputFile = outputFile;
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
        this.nonspecificEnrichSampler = new NonSpecificEnrichmentSampler(this.nonspecificEnrichParam1, this.nonspecificEnrichParam2);
        this.tretMethylationLevelSampler = new MethylationLevelSampler(this.methylationLevelParam1, this.methylationLevelParam2);
        this.ctrlMethylationLevelSampler = new MethylationLevelSampler(this.methylationLevelParam1, this.methylationLevelParam2);
        this.tretIPOverdispersionSampler = new OverdispersionSampler(this.ipOverdispersionShape, this.ipOverdispersionScale);
        this.tretINPUTOverdispersionSampler = new OverdispersionSampler(this.inputOverdispersionShape, this.inputOverdispersionScale);
        this.ctrlIPOverdispersionSampler = new OverdispersionSampler(this.ipOverdispersionShape, this.ipOverdispersionScale);
        this.ctrlINPUTOverdispersionSampler = new OverdispersionSampler(this.inputOverdispersionShape, this.inputOverdispersionScale);

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

        // non-peak region background expression of treatment group
        be = new BackgroundExpression(this.treatmentIPNonPeak, this.treatmentINPUTNonPeak);
        backgroundExpressionMean = be.geneBackgroundExp();
        backgroundExpressionStd = be.geneExpressionStd();
        this.tretNonPeakExpressionMean = new double[this.geneNumber];
        this.treatmentNonPeakExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            // shape and scale parameters of log-normal distribution
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = backgroundExpressionStd[geneIdx];
            this.tretNonPeakExpressionMean[geneIdx] = backgroundExpressionMean[geneIdx];
            this.treatmentNonPeakExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }
        this.tretIPNonPeakSizeFactors = be.getGlobalIPSizeFactor();
        this.tretINPUTNonPeakSizeFactors = be.getGlobalINPUTSizeFactor();

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

        // non-peak region background expression of control group
        be = new BackgroundExpression(this.controlIPNonPeak, this.controlINPUTNonPeak);
        backgroundExpressionMean = be.geneBackgroundExp();
        backgroundExpressionStd = be.geneExpressionStd();
        this.ctrlNonPeakExpressionMean = new double[this.geneNumber];
        this.controlNonPeakExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            // shape and scale parameters of log-normal distribution
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = backgroundExpressionStd[geneIdx];
            this.ctrlNonPeakExpressionMean[geneIdx] = backgroundExpressionMean[geneIdx];
            this.controlNonPeakExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }
        this.ctrlIPNonPeakSizeFactors = be.getGlobalIPSizeFactor();
        this.ctrlINPUTNonPeakSizeFactors = be.getGlobalINPUTSizeFactor();

        // treatment group methylation level, IP overdispersion, INPUT overdispersion of each gene, shape samplingTime × geneNumber
        double[] initMethLevel = new double[this.geneNumber];
        for (int i=0; i<this.geneNumber; i++) {
            initMethLevel[i] = 0.1;
        }

        double[] initNonspecific = new double[this.geneNumber];
        for (int i=0; i<this.geneNumber; i++) {
            initNonspecific[i] = 0.1;
        }
        ReadsExpectation re = new ReadsExpectation(this.treatmentIPReads, this.treatmentINPUTReads, this.treatmentIPNonPeak, this.treatmentINPUTNonPeak, initMethLevel, initNonspecific);
        this.tretIPReadsExpectation = re.getIPReadsExpectation();
        this.tretINPUTReadsExpectation = re.getINPUTReadsExpectation();
        this.tretIPNonPeakExpectation = re.getIPNonPeakExpectation();
        this.tretINPUTNonPeakExpectation = re.getINPUTNonPeakExpectation();

        re = new ReadsExpectation(this.controlIPReads, this.controlINPUTReads, this.controlIPNonPeak, this.controlINPUTNonPeak, initMethLevel, initNonspecific);
        this.ctrlIPReadsExpectation = re.getIPReadsExpectation();
        this.ctrlINPUTReadsExpectation = re.getINPUTReadsExpectation();
        this.ctrlIPNonPeakExpectation = re.getIPNonPeakExpectation();
        this.ctrlINPUTNonPeakExpectation = re.getINPUTNonPeakExpectation();
        initMethLevel = null;
        initNonspecific = null;
        re = null;
    }

    /**
     * pseudo prior distribution for IP and INPUT Overdispersion, methylation level, and background expression for each gene
     */
    private void selectModel() {
        this.bayesFactors = new ConcurrentHashMap<>();
        this.quantifyResult = new ConcurrentHashMap<>();
        ExecutorService threadPool = Executors.newFixedThreadPool(this.threadNum);
        CountDownLatch countDown = new CountDownLatch(this.geneNumber);
        CreateTask ct = (geneIdx -> () -> {
            try {
                 System.out.println(geneIdx+1 + "/" + this.geneNumber);
                int[][] geneReadsCounts;
                int[] tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                      tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak, selectResult;
                double[][] geneReadsExpectations;
                double[] tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                         tretIPNonPeakExpect, tretINPUTNonPeakExpect, ctrlIPNonPeakExpect, ctrlINPUTNonPeakExpect;
                SameMethylationLevelModel sameMethylationLevelModel;
                DiffMethylationLevelModel diffMethylationLevelModel;

                BackgroundExpressionSampler treatmentBkgExpSampler, treatmentNonPeakExpSampler, controlBkgExpSampler, controlNonPeakExpSampler;
                treatmentBkgExpSampler = this.treatmentBackgroundExpressionSamplers[geneIdx];
                controlBkgExpSampler = this.controlBackgroundExpressionSamplers[geneIdx];
                treatmentNonPeakExpSampler = this.treatmentNonPeakExpressionSamplers[geneIdx];
                controlNonPeakExpSampler = this.controlNonPeakExpressionSamplers[geneIdx];
                // data preparation
                geneReadsCounts = this.getGeneReads(geneIdx);
                tretIPReads = geneReadsCounts[0];
                tretINPUTReads = geneReadsCounts[1];
                ctrlIPReads = geneReadsCounts[2];
                ctrlINPUTReads = geneReadsCounts[3];
                tretIPNonPeak = geneReadsCounts[4];
                tretINPUTNonPeak = geneReadsCounts[5];
                ctrlIPNonPeak = geneReadsCounts[6];
                ctrlINPUTNonPeak = geneReadsCounts[7];

                geneReadsExpectations = this.getGeneExpectation(geneIdx);
                tretIPExpectation = geneReadsExpectations[0];
                tretINPUTExpectation = geneReadsExpectations[1];
                ctrlIPExpectation = geneReadsExpectations[2];
                ctrlINPUTExpectation = geneReadsExpectations[3];
                tretIPNonPeakExpect = geneReadsExpectations[4];
                tretINPUTNonPeakExpect = geneReadsExpectations[5];
                ctrlIPNonPeakExpect = geneReadsExpectations[6];
                ctrlINPUTNonPeakExpect = geneReadsExpectations[7];

                // initialize two models for a gene
                sameMethylationLevelModel = new SameMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                                          this.nonspecificEnrichSampler,
                                                                          this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                                          this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                                          treatmentBkgExpSampler, treatmentNonPeakExpSampler, controlBkgExpSampler, controlNonPeakExpSampler);
                diffMethylationLevelModel = new DiffMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                                          this.nonspecificEnrichSampler,
                                                                          this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                                          this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                                          treatmentBkgExpSampler, treatmentNonPeakExpSampler, controlBkgExpSampler, controlNonPeakExpSampler);
                // set gene reads count for each model
                sameMethylationLevelModel.setTreatmentControlGeneReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                                                                       tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak,
                                                                       tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                       tretIPNonPeakExpect, tretINPUTNonPeakExpect, ctrlIPNonPeakExpect, ctrlINPUTNonPeakExpect);
                diffMethylationLevelModel.setTreatmentControlGeneReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                                                                       tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak,
                                                                       tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                       tretIPNonPeakExpect, tretINPUTNonPeakExpect, ctrlIPNonPeakExpect, ctrlINPUTNonPeakExpect);

                // set individual size factors
                sameMethylationLevelModel.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors,
                                                         this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
                diffMethylationLevelModel.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors,
                                                         this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
                // first time sampling, set sampling list for each model
                this.setModelSamplingList(sameMethylationLevelModel, diffMethylationLevelModel, geneIdx);
                selectResult = this.sampling(sameMethylationLevelModel, diffMethylationLevelModel, false);

                // estimate pseudo prior distribution parameters and start the second time sampling
                // sameMethylationLevelModel.setModelParameterPseudoPrior();
                // diffMethylationLevelModel.setModelParameterPseudoPrior();
                // this.setModelSamplingList(sameMethylationLevelModel, diffMethylationLevelModel, geneIdx);
                // selectResult = this.sampling(sameMethylationLevelModel, diffMethylationLevelModel, true);
                // System.out.println(Arrays.stream(selectResult).skip(this.burnIn).sum());

                // calculate same methylation level model posterior density and diff methylation level model posterior density
                // double diffMethModelProba = this.differentiation(selectResult);
                // double sameMethModelProba = 1 - diffMethModelProba;

                // ModelSelection finalModel = (diffMethModelProba-0.5<0.00001)? sameMethylationLevelModel: diffMethylationLevelModel;
                // double[] finalModelQuantify = this.quantify(finalModel);
                double[] model0Quantify = this.quantify(sameMethylationLevelModel);
                double[] model1Quantify = this.quantify(diffMethylationLevelModel);

                // this.quantifyResult.put(geneIdx, finalModelQuantify);
                // calculate Bayes factor
                // double bayesFactor = this.calcLogBayesFactorLaplaceApproximation(model0Quantify, model1Quantify, sameMethylationLevelModel, diffMethylationLevelModel, geneIdx);
                double bayesFactor = this.calcBayesFactorLaplaceMetropolis(sameMethylationLevelModel, diffMethylationLevelModel);
                this.bayesFactors.put(geneIdx, bayesFactor);
                if (bayesFactor - 3 > 0.00001)
                    this.quantifyResult.put(geneIdx, model1Quantify);
                else
                    this.quantifyResult.put(geneIdx, model0Quantify);

                // help GC
                sameMethylationLevelModel.helpGC();
                diffMethylationLevelModel.helpGC();
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                countDown.countDown();
            }
        });
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx ++) {
            Runnable task = ct.getTask(geneIdx);
            threadPool.submit(task);
        }

        try {
            countDown.await();
            threadPool.shutdown();
            if (!threadPool.awaitTermination(1000, TimeUnit.MILLISECONDS))
                threadPool.shutdownNow();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
        } finally {
            threadPool.shutdownNow();
        }
    }

    /**
     * sampling process
     */
    private int[] sampling(SameMethylationLevelModel sameMethylationLevelModel, DiffMethylationLevelModel diffMethylationLevelModel, boolean usePseudoPrior) {
        int modelType;
        double sameMethModelLogProba, diffMethModelLogProba,
               sameMethFullConditionProba, diffMethFullConditionProba;
        // model type initialization, 0 denotes as same methylation model, 1 denotes as different methylation model
        int[] modelSelectionResult = new int[this.samplingTime];
        modelSelectionResult[0] = 0;

        for (int t=1; t<this.samplingTime; t++) {
            // log-likelihood + log-prior
            sameMethFullConditionProba = sameMethylationLevelModel.paramsLogFullConditionalProbability();
            diffMethFullConditionProba = diffMethylationLevelModel.paramsLogFullConditionalProbability();

            modelType = modelSelectionResult[t-1];
            sameMethFullConditionProba = sameMethylationLevelModel.iterate(t, sameMethFullConditionProba,true, usePseudoPrior);
            diffMethFullConditionProba = diffMethylationLevelModel.iterate(t, diffMethFullConditionProba,true, usePseudoPrior);

            // model full-conditional posterior probability = log-likelihood + log-prior + log-pseudo prior + log-model prior
            if (modelType == 0) { // current model is same methylation model
                if (t >= this.burnIn)
                    sameMethylationLevelModel.recordCurrentIterationParams();
                sameMethModelLogProba = sameMethFullConditionProba + diffMethylationLevelModel.paramsPseudoDistributionProbability(usePseudoPrior) + Math.log(this.sameMethLevelModelPriorProba);
                diffMethFullConditionProba = diffMethylationLevelModel.paramsLogFullConditionalProbability();
                diffMethModelLogProba = diffMethFullConditionProba + sameMethylationLevelModel.paramsPseudoDistributionProbability(usePseudoPrior) + Math.log(this.diffMethLevelModelPriorProba);
            } else {    // current model is different methylation model
                if (t >= this.burnIn)
                    diffMethylationLevelModel.recordCurrentIterationParams();
                sameMethFullConditionProba = sameMethylationLevelModel.paramsLogFullConditionalProbability();
                sameMethModelLogProba = sameMethFullConditionProba + diffMethylationLevelModel.paramsPseudoDistributionProbability(usePseudoPrior) + Math.log(this.sameMethLevelModelPriorProba);
                diffMethModelLogProba = diffMethFullConditionProba + sameMethylationLevelModel.paramsPseudoDistributionProbability(usePseudoPrior) + Math.log(this.diffMethLevelModelPriorProba);
            }
            modelSelectionResult[t] = this.modelSelector.getSamplingRes(sameMethModelLogProba, diffMethModelLogProba);
        }

        return modelSelectionResult;
    }

    /**
     * treatment and control group methylation, IP overdispersion, INPUT overdispersion, background expression quantification
     * use quantification model Quantification.MethylationQuantification
     * @param model model need to be quantified
     * @return double[] {treatment methylation, control methylation, treatment IP overdispersion, treatment INPUT overdispersion,
     *                   control IP overdispersion, control INPUT overdispersion, treatment background expression, control background expression}
     */
    private double[] quantify(ModelSelection model) {
        int remainTime = this.samplingTime - this.burnIn;
        double[] treatmentMethRemainValue = Arrays.stream(model.treatmentMethylationLevel).skip(this.burnIn).sorted().toArray();
        double[] controlMethRemainValue = Arrays.stream(model.controlMethylationLevel).skip(this.burnIn).sorted().toArray();
        double[] nonspecificEnrichRemainValue = Arrays.stream(model.nonspecificEnrichment).skip(this.burnIn).sorted().toArray();
        double[] treatmentIPOverdispersionRemainValue = Arrays.stream(model.treatmentIPOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] treatmentINPUTOverdispersionRemainValue = Arrays.stream(model.treatmentINPUTOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] controlIPOverdispersionRemainValue = Arrays.stream(model.controlIPOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] controlINPUTOverdispersionRemainValue = Arrays.stream(model.controlINPUTOverdispersion).skip(this.burnIn).sorted().toArray();
//        double[] treatmentBkgExpRemainValue = Arrays.stream(model.treatmentBkgExp).skip(this.burnIn).sorted().toArray();
//        double[] controlBkgExpRemainValue = Arrays.stream(model.controlBkgExp).skip(this.burnIn).sorted().toArray();
//        double[] treatmentNonPeakExpRemainValue = Arrays.stream(model.treatmentNonPeakBkgExp).skip(this.burnIn).sorted().toArray();
//        double[] controlNonPeakExpRemainValue = Arrays.stream(model.controlNonPeakBkgExp).skip(this.burnIn).sorted().toArray();

        int medianIdx = remainTime / 2;
        double tretMeth, ctrlMeth, nonspecificEnrich, tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion,
               ctrlINPUTOverdispersion; // , tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp
        if (remainTime%2==0) {
            tretMeth = (treatmentMethRemainValue[medianIdx] + treatmentMethRemainValue[medianIdx+1]) / 2;
            ctrlMeth = (controlMethRemainValue[medianIdx] + controlMethRemainValue[medianIdx+1]) / 2;
            nonspecificEnrich = (nonspecificEnrichRemainValue[medianIdx] + nonspecificEnrichRemainValue[medianIdx+1]) / 2;
            tretIPOverdispersion = (treatmentIPOverdispersionRemainValue[medianIdx] + treatmentIPOverdispersionRemainValue[medianIdx+1]) / 2;
            tretINPUTOverdispersion = (treatmentINPUTOverdispersionRemainValue[medianIdx] + treatmentINPUTOverdispersionRemainValue[medianIdx+1]) / 2;
            ctrlIPOverdispersion = (controlIPOverdispersionRemainValue[medianIdx] + controlIPOverdispersionRemainValue[medianIdx+1]) / 2;
            ctrlINPUTOverdispersion = (controlINPUTOverdispersionRemainValue[medianIdx] + controlINPUTOverdispersionRemainValue[medianIdx+1]) / 2;
//            tretBkgExp = (treatmentBkgExpRemainValue[medianIdx] + treatmentBkgExpRemainValue[medianIdx+1]) / 2;
//            ctrlBkgExp = (controlBkgExpRemainValue[medianIdx] + controlBkgExpRemainValue[medianIdx+1]) / 2;
//            tretNonPeakExp = (treatmentNonPeakExpRemainValue[medianIdx] + treatmentNonPeakExpRemainValue[medianIdx+1]) / 2;
//            ctrlNonPeakExp = (controlNonPeakExpRemainValue[medianIdx] + controlNonPeakExpRemainValue[medianIdx+1]) / 2;
        } else {
            tretMeth = treatmentMethRemainValue[medianIdx];
            ctrlMeth = controlMethRemainValue[medianIdx];
            nonspecificEnrich = nonspecificEnrichRemainValue[medianIdx];
            tretIPOverdispersion = treatmentIPOverdispersionRemainValue[medianIdx];
            tretINPUTOverdispersion = treatmentINPUTOverdispersionRemainValue[medianIdx];
            ctrlIPOverdispersion = controlIPOverdispersionRemainValue[medianIdx];
            ctrlINPUTOverdispersion = controlINPUTOverdispersionRemainValue[medianIdx];
//            tretBkgExp = treatmentBkgExpRemainValue[medianIdx];
//            ctrlBkgExp = controlBkgExpRemainValue[medianIdx];
//            tretNonPeakExp = treatmentNonPeakExpRemainValue[medianIdx];
//            ctrlNonPeakExp = controlNonPeakExpRemainValue[medianIdx];
        }

        return new double[] {tretMeth, ctrlMeth, nonspecificEnrich,
                             tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion};
                            // , tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp
    }

    /**
     * Deprecated
     * probability of variant methylation level between treatment and control
     * @return probability
     */
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
    private double calcLogBayesFactorLaplaceApproximation(double[] sameMethModelQuantify, double[] diffMethModelQuantify,
                                                          ModelSelection sameMethModel, ModelSelection diffMethModel, int geneIdx) {
        int[] tretIPReads = sameMethModel.treatmentIPReads;
        int[] tretINPUTReads = sameMethModel.treatmentINPUTReads;
        int[] ctrlIPReads = sameMethModel.controlIPReads;
        int[] ctrlINPUTReads = sameMethModel.controlINPUTReads;
        int[] tretIPNonPeak = sameMethModel.treatmentIPNonPeakReads;
        int[] tretINPUTNonPeak = sameMethModel.treatmentINPUTNonPeakReads;
        int[] ctrlIPNonPeak = sameMethModel.controlIPNonPeakReads;
        int[] ctrlINPUTNonPeak = sameMethModel.controlINPUTNonPeakReads;
        double tretPeakBkgExp = this.tretGeneBackgroundExpression[geneIdx];
        double tretNonPeakBkgExp= this.tretNonPeakExpressionMean[geneIdx];
        double ctrlPeakBkgExp = this.ctrlGeneBackgroundExpression[geneIdx];
        double ctrlNonPeakBkgExp= this.ctrlNonPeakExpressionMean[geneIdx];
        double model0LogLikelihood = this.modelLogLikelihood(sameMethModelQuantify, sameMethModel);
        double model1LofLikelihood = this.modelLogLikelihood(diffMethModelQuantify, diffMethModel);
        double model0LogPrior = this.modelLogPrior(sameMethModelQuantify, sameMethModel);
        double model1LogPrior = this.modelLogPrior(diffMethModelQuantify, diffMethModel);
        double[] tretSizeFactor = this.sampleSizeFactors(this.treatmentIPReads, this.treatmentINPUTReads);
        double[] ctrlSizeFactor = this.sampleSizeFactors(this.controlIPReads, this.controlINPUTReads);
        BayesFactorCalculation bfc = new BayesFactorCalculation(model0LogLikelihood, model1LofLikelihood, model0LogPrior, model1LogPrior,
                                                                sameMethModelQuantify, diffMethModelQuantify, tretPeakBkgExp, ctrlPeakBkgExp, tretNonPeakBkgExp, ctrlNonPeakBkgExp,
                                                                this.ipOverdispersionShape, this.methylationLevelParam1, this.methylationLevelParam2);
        bfc.setReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                     tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak);
        bfc.setSizeFactors(tretSizeFactor, ctrlSizeFactor);
        return bfc.calcBayesFactor();
    }

    /**
     * calculate Bayes factor with Laplace-Metropolis estimator
     * @return Bayes factor
     */
    private double calcBayesFactorLaplaceMetropolis(ModelSelection sameMethModel, ModelSelection diffMethModel) {
        double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                 treatmentMethylationLevel, controlMethylationLevel, nonspecificEnrichment;

        treatmentIPOverdispersion = sameMethModel.treatmentIPOverdispersion;
        treatmentINPUTOverdispersion = sameMethModel.treatmentINPUTOverdispersion;
        controlIPOverdispersion = sameMethModel.controlIPOverdispersion;
        controlINPUTOverdispersion = sameMethModel.controlINPUTOverdispersion;
        treatmentMethylationLevel = sameMethModel.treatmentMethylationLevel;
        nonspecificEnrichment = sameMethModel.nonspecificEnrichment;
        CovarianceMatrix model0Covariance = new CovarianceMatrix(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                                 controlIPOverdispersion, controlINPUTOverdispersion,
                                                                 treatmentMethylationLevel, null,
                                                                 nonspecificEnrichment, this.burnIn);
        double[][] model0CovarianceMatrix = model0Covariance.getCovarianceMatrix();
        RealMatrix model0CovMatrix = new Array2DRowRealMatrix(model0CovarianceMatrix);
        double model0Determinate = new LUDecomposition(model0CovMatrix).getDeterminant();
        BestModelParams model0BestParams = new BestModelParams(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                               controlIPOverdispersion, controlINPUTOverdispersion,
                                                               treatmentMethylationLevel, null,
                                                               nonspecificEnrichment, sameMethModel, this.burnIn);
        model0BestParams.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors,
                                        this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
        double model0MaximumPosterior = model0BestParams.getMaximumPosterior();

        treatmentIPOverdispersion = diffMethModel.treatmentIPOverdispersion;
        treatmentINPUTOverdispersion = diffMethModel.treatmentINPUTOverdispersion;
        controlIPOverdispersion = diffMethModel.controlIPOverdispersion;
        controlINPUTOverdispersion = diffMethModel.controlINPUTOverdispersion;
        treatmentMethylationLevel = diffMethModel.treatmentMethylationLevel;
        controlMethylationLevel = diffMethModel.controlMethylationLevel;
        nonspecificEnrichment = diffMethModel.nonspecificEnrichment;
        CovarianceMatrix model1Covariance = new CovarianceMatrix(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                                 controlIPOverdispersion, controlINPUTOverdispersion,
                                                                 treatmentMethylationLevel, controlMethylationLevel,
                                                                 nonspecificEnrichment, this.burnIn);
        double[][] model1CovarianceMatrix = model1Covariance.getCovarianceMatrix();
        RealMatrix model1CovMatrix = new Array2DRowRealMatrix(model1CovarianceMatrix);
        double model1Determinate = new LUDecomposition(model1CovMatrix).getDeterminant();
        BestModelParams model1BestParams = new BestModelParams(treatmentIPOverdispersion, treatmentINPUTOverdispersion,
                                                               controlIPOverdispersion, controlINPUTOverdispersion,
                                                               treatmentMethylationLevel, controlMethylationLevel,
                                                               nonspecificEnrichment, sameMethModel, this.burnIn);
        model1BestParams.setSizeFactors(this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors,
                                        this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
        double model1MaximumPosterior = model1BestParams.getMaximumPosterior();

        int model0ParamNum = 6, model1ParamNum = 7;

        // p(D|M=0) and p(D|M=1)
        double model0Proba = Math.pow(2*Math.PI, model0ParamNum*0.5) * Math.pow(model0Determinate, 0.5) * model0MaximumPosterior;
        double model1Proba = Math.pow(2*Math.PI, model1ParamNum*0.5) * Math.pow(model1Determinate, 0.5) * model1MaximumPosterior;

        return model1Proba / model0Proba;
    }

    /**
     * size factors of each individual
     * @param ipReads IP reads, shape individualNumber × geneNumber
     * @param inputReads INPUT reads, shape individualNumber × geneNumber
     * @return samples size factors
     */
    private double[] sampleSizeFactors(int[][] ipReads, int[][] inputReads) {
        BackgroundExpression be = new BackgroundExpression(ipReads, inputReads);
        return be.getGlobalIPSizeFactor();
        // only for simulation data test
        // return new double[] {1, 1, 1};
    }

    /**
     * calculate log-likelihood with quantification result
     * @param quantify quantification result
     * @param model model
     * @return log-likelihood
     */
    private double modelLogLikelihood(double[] quantify, ModelSelection model) {
        double tretMeth = quantify[0];
        double ctrlMeth = quantify[1];
        double nonspecificEnrich = quantify[2];
        double tretIPOverdispersion = quantify[3];
        double tretINPUTOverdispersion = quantify[4];
        double ctrlIPOverdispersion = quantify[5];
        double ctrlINPUTOverdispersion = quantify[6];
        double[] newTretIPReadsExpectation, newTretIPNonPeakExpectation, newTretINPUTReadsExpectation, newTretINPUTNonPeakExpectation,
                 newCtrlIPReadsExpectation, newCtrlIPNonPeakExpectation, newCtrlINPUTReadsExpectation, newCtrlINPUTNonPeakExpectation;
        double[][] newExpectations;

        newExpectations = model.renewReadsExpectationViaNonspecificEnrich(true, new double[]{tretMeth}, new double[]{nonspecificEnrich},
                                                                          this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors);
        newTretIPReadsExpectation = newExpectations[0];
        newTretIPNonPeakExpectation = newExpectations[1];
        newTretINPUTReadsExpectation = newExpectations[2];
        newTretINPUTNonPeakExpectation = newExpectations[3];

        newExpectations = model.renewReadsExpectationViaNonspecificEnrich(false, new double[]{ctrlMeth}, new double[]{nonspecificEnrich},
                                                                           this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
        newCtrlIPReadsExpectation = newExpectations[0];
        newCtrlIPNonPeakExpectation = newExpectations[1];
        newCtrlINPUTReadsExpectation = newExpectations[2];
        newCtrlINPUTNonPeakExpectation = newExpectations[3];

        return model.logLikelihood(newTretIPReadsExpectation, newTretINPUTReadsExpectation, newCtrlIPReadsExpectation, newCtrlINPUTReadsExpectation,
                                   newTretIPNonPeakExpectation, newTretINPUTNonPeakExpectation, newCtrlIPNonPeakExpectation, newCtrlINPUTNonPeakExpectation,
                                   tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
    }

    /**
     * calculate log-prior with quantification result
     * @param quantify quantification result
     * @param model model
     * @return log-likelihood
     */
    private double modelLogPrior(double[] quantify, ModelSelection model) {
        double tretMeth = quantify[0];
        double ctrlMeth = quantify[1];
        double nonspecificEnrich = quantify[2];
        double tretIPOverdispersion = quantify[3];
        double tretINPUTOverdispersion = quantify[4];
        double ctrlIPOverdispersion = quantify[5];
        double ctrlINPUTOverdispersion = quantify[6];

        return model.logPriority(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                 tretMeth, ctrlMeth, nonspecificEnrich, 0,0,0,0);
    }

    /**
     * treatment and control group IP and INPUT reads count
     * @param geneIdx gene index
     * @return reads count, shape 1 × individualNumber
     */
    private int[][] getGeneReads(int geneIdx) {
        int[] tretIPReads = new int[this.tretIndividualNumber], tretINPUTReads = new int[this.tretIndividualNumber],
              ctrlIPReads = new int[this.ctrlIndividualNumber], ctrlINPUTReads = new int[this.ctrlIndividualNumber],
              tretIPNonPeak = new int[this.tretIndividualNumber], tretINPUTNonPeak = new int[this.tretIndividualNumber],
              ctrlIPNonPeak = new int[this.ctrlIndividualNumber], ctrlINPUTNonPeak = new int[this.ctrlIndividualNumber];

        for (int individualIdx=0; individualIdx<this.tretIndividualNumber; individualIdx++) {
            tretIPReads[individualIdx] = this.treatmentIPReads[individualIdx][geneIdx];
            tretINPUTReads[individualIdx] = this.treatmentINPUTReads[individualIdx][geneIdx];
            tretIPNonPeak[individualIdx] = this.treatmentIPNonPeak[individualIdx][geneIdx];
            tretINPUTNonPeak[individualIdx] = this.treatmentINPUTNonPeak[individualIdx][geneIdx];
        }

        for (int individualIdx=0; individualIdx<this.ctrlIndividualNumber; individualIdx++) {
            ctrlIPReads[individualIdx] = this.controlIPReads[individualIdx][geneIdx];
            ctrlINPUTReads[individualIdx] = this.controlINPUTReads[individualIdx][geneIdx];
            ctrlIPNonPeak[individualIdx] = this.controlIPNonPeak[individualIdx][geneIdx];
            ctrlINPUTNonPeak[individualIdx] = this.controlINPUTNonPeak[individualIdx][geneIdx];
        }

        return new int[][] {tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                            tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak};
    }

    /**
     * treatment and control group IP and INPUT reads count expectations
     * @param geneIdx gene index
     * @return reads count expectations, shape 1 × individualNumber
     */
    private double[][] getGeneExpectation(int geneIdx) {
        double[] tretIPReads = new double[this.tretIndividualNumber], tretINPUTReads = new double[this.tretIndividualNumber],
                 ctrlIPReads = new double[this.ctrlIndividualNumber], ctrlINPUTReads = new double[this.ctrlIndividualNumber],
                 tretIPNonPeak = new double[this.tretIndividualNumber], tretINPUTNonPeak = new double[this.tretIndividualNumber],
                 ctrlIPNonPeak = new double[this.ctrlIndividualNumber], ctrlINPUTNonPeak = new double[this.ctrlIndividualNumber];

        for (int individualIdx=0; individualIdx<this.tretIndividualNumber; individualIdx++) {
            tretIPReads[individualIdx] = this.tretIPReadsExpectation[individualIdx][geneIdx];
            tretINPUTReads[individualIdx] = this.tretINPUTReadsExpectation[individualIdx][geneIdx];
            tretIPNonPeak[individualIdx] = this.tretIPNonPeakExpectation[individualIdx][geneIdx];
            tretINPUTNonPeak[individualIdx] = this.tretINPUTNonPeakExpectation[individualIdx][geneIdx];
        }

        for (int individualIdx=0; individualIdx<this.ctrlIndividualNumber; individualIdx++) {
            ctrlIPReads[individualIdx] = this.ctrlIPReadsExpectation[individualIdx][geneIdx];
            ctrlINPUTReads[individualIdx] = this.ctrlINPUTReadsExpectation[individualIdx][geneIdx];
            ctrlIPNonPeak[individualIdx] = this.ctrlIPNonPeakExpectation[individualIdx][geneIdx];
            ctrlINPUTNonPeak[individualIdx] = this.ctrlINPUTNonPeakExpectation[individualIdx][geneIdx];
        }

        return new double[][] {tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                               tretIPNonPeak, tretINPUTNonPeak, ctrlIPNonPeak, ctrlINPUTNonPeak};
    }

    /**
     * set sampling list for each model
     * @param geneIdx gene index
     */
    private void setModelSamplingList(SameMethylationLevelModel sameMethylationLevelModel, DiffMethylationLevelModel diffMethylationLevelModel, int geneIdx) {
        // here m1, m2 represent same methylation level model and variant methylation level model, respectively
        double[] m1TretMethLevel, m2TretMethLevel, m1TretBkgExp, m2TretBkgExp, m1TretIPOverdispersion, m1TretINPUTOverdispersion, m2TretIPOverdispersion, m2TretINPUTOverdispersion,
                m1CtrlMethLevel, m2CtrlMethLevel, m1CtrlBkgExp, m2CtrlBkgExp, m1CtrlIPOverdispersion, m1CtrlINPUTOverdispersion, m2CtrlIPOverdispersion, m2CtrlINPUTOverdispersion,
                m1TretNonPeakExp, m1CtrlNonPeakExp, m2TretNonPeakExp, m2CtrlNonPeakExp, m1NonspecificEnrich, m2NonspecificEnrich;
        // set sampling list for each model
        m1TretMethLevel = new double[this.samplingTime];
        m2TretMethLevel = new double[this.samplingTime];
        m1TretMethLevel[0] = 0.01;
        m2TretMethLevel[0] = 0.01;
        m1CtrlMethLevel = new double[this.samplingTime];
        m2CtrlMethLevel = new double[this.samplingTime];
        m1CtrlMethLevel[0] = 0.01;
        m2CtrlMethLevel[0] = 0.01;

        m1TretBkgExp = new double[this.samplingTime];
        m2TretBkgExp = new double[this.samplingTime];
        m1TretBkgExp[0] = this.tretGeneBackgroundExpression[geneIdx];
        m2TretBkgExp[0] = this.tretGeneBackgroundExpression[geneIdx];
        m1CtrlBkgExp = new double[this.samplingTime];
        m2CtrlBkgExp = new double[this.samplingTime];
        m1CtrlBkgExp[0] = this.ctrlGeneBackgroundExpression[geneIdx];
        m2CtrlBkgExp[0] = this.ctrlGeneBackgroundExpression[geneIdx];

        m1TretNonPeakExp = new double[this.samplingTime];
        m2TretNonPeakExp = new double[this.samplingTime];
        m1TretNonPeakExp[0] = this.tretNonPeakExpressionMean[geneIdx];
        m2TretNonPeakExp[0] = this.tretNonPeakExpressionMean[geneIdx];
        m1CtrlNonPeakExp = new double[this.samplingTime];
        m2CtrlNonPeakExp = new double[this.samplingTime];
        m1CtrlNonPeakExp[0] = this.ctrlNonPeakExpressionMean[geneIdx];
        m2CtrlNonPeakExp[0] = this.ctrlNonPeakExpressionMean[geneIdx];

        m1TretIPOverdispersion = new double[this.samplingTime];
        m2TretIPOverdispersion = new double[this.samplingTime];
        m1TretIPOverdispersion[0] = 0.1;
        m2TretIPOverdispersion[0] = 0.1;
        m1TretINPUTOverdispersion = new double[this.samplingTime];
        m2TretINPUTOverdispersion = new double[this.samplingTime];
        m1TretINPUTOverdispersion[0] = 0.1;
        m2TretINPUTOverdispersion[0] = 0.1;

        m1CtrlIPOverdispersion = new double[this.samplingTime];
        m2CtrlIPOverdispersion = new double[this.samplingTime];
        m1CtrlIPOverdispersion[0] = 0.1;
        m2CtrlIPOverdispersion[0] = 0.1;
        m1CtrlINPUTOverdispersion = new double[this.samplingTime];
        m2CtrlINPUTOverdispersion = new double[this.samplingTime];
        m1CtrlINPUTOverdispersion[0] = 0.1;
        m2CtrlINPUTOverdispersion[0] = 0.1;

        m1NonspecificEnrich = new double[this.samplingTime];
        m2NonspecificEnrich = new double[this.samplingTime];
        m1NonspecificEnrich[0] = 0.1;
        m2NonspecificEnrich[0] = 0.1;

        sameMethylationLevelModel.setSamplingList(m1TretIPOverdispersion, m1TretINPUTOverdispersion,
                                                  m1CtrlIPOverdispersion, m1CtrlINPUTOverdispersion,
                                                  m1TretBkgExp, m1CtrlBkgExp, m1TretNonPeakExp, m1CtrlNonPeakExp,
                                                  m1TretMethLevel, m1CtrlMethLevel, m1NonspecificEnrich);
        diffMethylationLevelModel.setSamplingList(m2TretIPOverdispersion, m2TretINPUTOverdispersion,
                                                  m2CtrlIPOverdispersion, m2CtrlINPUTOverdispersion,
                                                  m2TretBkgExp, m2CtrlBkgExp, m2TretNonPeakExp, m2CtrlNonPeakExp,
                                                  m2TretMethLevel, m2CtrlMethLevel, m2NonspecificEnrich);
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
            bfw.write("#geneIdx\tBF\ttretMeth\tctrlMeth\tnonspecificEnrich\ttretIPOverdispersion\ttretINPUTOverdispersion\tctrlIPOverdispersion\tctrlINPUTOverdispersion");//\ttretBkgExp\tctrlBkgExp\ttretNonPeakExp\tctrlNonPeakExp
            bfw.newLine();
            String line, records, bf;
            double[] quantify;
            for (Integer geneIdx: sortedList) {
                quantify = this.quantifyResult.get(geneIdx);
                records = Arrays.stream(quantify).mapToObj(x -> (this.df.format((Double) x))).collect(Collectors.joining("\t"));
                bf = this.df.format(this.bayesFactors.get(geneIdx));
                line = geneIdx+1 + "\t" + bf + "\t" + records;
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
