package DifferentialMethylation;

import Quantification.BackgroundExpressionSampler;
import Quantification.MethylationLevelSampler;
import Quantification.OverdispersionSampler;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.ReadsExpectation;

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
    private int individualNumber, geneNumber, samplingTime, burnIn, threadNum;
    private int[][] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads;
    // distribution parameters
    private double ipOverdispersionShape, ipOverdispersionScale, inputOverdispersionShape, inputOverdispersionScale,
                   methylationLevelParam1, methylationLevelParam2, sameMethLevelModelPriorProba, diffMethLevelModelPriorProba;
    private double[][] tretIPReadsExpectation, tretINPUTReadsExpectation, ctrlIPReadsExpectation, ctrlINPUTReadsExpectation;
    // MH sampling components, shape geneNumber × samplingTime
    private double[] tretGeneBackgroundExpression, ctrlGeneBackgroundExpression;
    private String outputFile;
    private OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler;
    private MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler;
    private BackgroundExpressionSampler[] treatmentBackgroundExpressionSamplers, controlBackgroundExpressionSamplers;
    private ModelSelector modelSelector;
    private ConcurrentHashMap<Integer, Double> bayesFactors;
    private ConcurrentHashMap<Integer, double[]> quantifyResult, modelProbabilities;
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
    public MethylationDifference(int[][] treatmentIPReads, int[][] treatmentINPUTReads,
                                 int[][] controlIPReads, int[][] controlINPUTReads,
                                 double ipOverdispersionShape, double ipOverdispersionScale,
                                 double inputOverdispersionShape, double inputOverdispersionScale,
                                 double methylationLevelParam1, double methylationLevelParam2,
                                 double sameMethLevelModelPriorProba, double diffMethLevelModelPriorProba,
                                 int samplingTime, int burnIn, String outputFile, int threadNumber) {
        // checking input parameters
        assert Math.abs(sameMethLevelModelPriorProba+diffMethLevelModelPriorProba-1) < 0.00001;
        assert treatmentINPUTReads.length == treatmentIPReads.length && controlINPUTReads.length == controlIPReads.length;
        assert treatmentINPUTReads.length == controlINPUTReads.length;
        assert treatmentINPUTReads[0].length == treatmentIPReads[0].length && controlINPUTReads[0].length == controlIPReads[0].length;
        assert treatmentINPUTReads[0].length == controlINPUTReads[0].length;
        assert samplingTime > burnIn;

        this.controlINPUTReads = controlINPUTReads;
        this.controlIPReads = controlIPReads;
        this.treatmentINPUTReads = treatmentINPUTReads;
        this.treatmentIPReads = treatmentIPReads;
        this.individualNumber = treatmentINPUTReads.length;
        this.geneNumber = treatmentINPUTReads[0].length;
        this.ipOverdispersionShape = ipOverdispersionShape;
        this.ipOverdispersionScale = ipOverdispersionScale;
        this.inputOverdispersionShape = inputOverdispersionShape;
        this.inputOverdispersionScale = inputOverdispersionScale;
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

        // background expression of treatment group
        BackgroundExpression be = new BackgroundExpression(this.treatmentIPReads, this.treatmentINPUTReads);
        double[] backgroundExpressionMean = be.geneBackgroundExp(); // background expression expectation, shape 1 × geneNumber
        double[] backgroundExpressionStd = be.geneExpressionStd(); // background expression standard deviation, shape 1 × geneNumber
        double scale, shape;
        this.tretGeneBackgroundExpression = new double[this.geneNumber];
        this.treatmentBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            // shape and scale parameters of log-normal distribution
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = Math.log(backgroundExpressionStd[geneIdx]);
            this.tretGeneBackgroundExpression[geneIdx] = backgroundExpressionMean[geneIdx];
            this.treatmentBackgroundExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }

        // background expression of control group
        be = new BackgroundExpression(this.controlIPReads, this.controlINPUTReads);
        backgroundExpressionMean = be.geneBackgroundExp();
        backgroundExpressionStd = be.geneExpressionStd();
        this.ctrlGeneBackgroundExpression = new double[this.geneNumber];
        this.controlBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            scale = Math.log(backgroundExpressionMean[geneIdx]);
            shape = Math.log(backgroundExpressionStd[geneIdx]);
            this.ctrlGeneBackgroundExpression[geneIdx] = backgroundExpressionMean[geneIdx];
            this.controlBackgroundExpressionSamplers[geneIdx] = new BackgroundExpressionSampler(scale, shape);
        }

        // treatment group methylation level, IP overdispersion, INPUT overdispersion of each gene, shape samplingTime × geneNumber
        double[] initMethLevel = new double[this.geneNumber];
        for (int i=0; i<this.geneNumber; i++) {
            initMethLevel[i] = 0.1;
        }
        ReadsExpectation re = new ReadsExpectation(this.treatmentIPReads, this.treatmentINPUTReads, initMethLevel);
        this.tretIPReadsExpectation = re.getIPReadsExepectation();
        this.tretINPUTReadsExpectation = re.getINPUTReadsExpectation();

        re = new ReadsExpectation(this.controlIPReads, this.controlINPUTReads, initMethLevel);
        this.ctrlIPReadsExpectation = re.getIPReadsExepectation();
        this.ctrlINPUTReadsExpectation = re.getINPUTReadsExpectation();
        initMethLevel = null;
        re = null;
    }

    /**
     * pseudo prior distribution for IP and INPUT Overdispersion, methylation level, and background expression for each gene
     */
    private void selectModel() {
        this.bayesFactors = new ConcurrentHashMap<>();
        this.quantifyResult = new ConcurrentHashMap<>();
        this.modelProbabilities = new ConcurrentHashMap<>();
        ExecutorService threadPool = Executors.newFixedThreadPool(this.threadNum);
        CountDownLatch countDown = new CountDownLatch(this.geneNumber);
        CreateTask ct = (geneIdx -> {
            return () -> {
                System.out.println(geneIdx);
                int[][] geneReadsCounts;
                int[] tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads, selectResult;
                double[][] geneReadsExpectations;
                double[] tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation;
                SameMethylationLevelModel sameMethylationLevelModel;
                DiffMethylationLevelModel diffMethylationLevelModel;

                BackgroundExpressionSampler treatmentBkgExpSampler, controlBkgExpSampler;
                treatmentBkgExpSampler = this.treatmentBackgroundExpressionSamplers[geneIdx];
                controlBkgExpSampler = this.controlBackgroundExpressionSamplers[geneIdx];
                // data preparation
                geneReadsCounts = this.getGeneReads(geneIdx);
                tretIPReads = geneReadsCounts[0];
                tretINPUTReads = geneReadsCounts[1];
                ctrlIPReads = geneReadsCounts[2];
                ctrlINPUTReads = geneReadsCounts[3];
                geneReadsCounts = null;
                geneReadsExpectations = this.getGeneExpectation(geneIdx);
                tretIPExpectation = geneReadsExpectations[0];
                tretINPUTExpectation = geneReadsExpectations[1];
                ctrlIPExpectation = geneReadsExpectations[2];
                ctrlINPUTExpectation = geneReadsExpectations[3];
                geneReadsExpectations = null;
                // initialize two models for a gene
                sameMethylationLevelModel = new SameMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                                          this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                                          this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                                          treatmentBkgExpSampler, controlBkgExpSampler);
                diffMethylationLevelModel = new DiffMethylationLevelModel(this.tretMethylationLevelSampler, this.ctrlMethylationLevelSampler,
                                                                          this.tretIPOverdispersionSampler, this.tretINPUTOverdispersionSampler,
                                                                          this.ctrlIPOverdispersionSampler, this.ctrlINPUTOverdispersionSampler,
                                                                          treatmentBkgExpSampler, controlBkgExpSampler);

                // set gene reads count for each model
                sameMethylationLevelModel.setTreatmentControlGeneReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                        tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation);
                diffMethylationLevelModel.setTreatmentControlGeneReads(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                        tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation);
                // first time sampling, set sampling list for each model
                this.setModelSamplingList(sameMethylationLevelModel, diffMethylationLevelModel, geneIdx);
                this.sampling(sameMethylationLevelModel, diffMethylationLevelModel, false);
                // estimate pseudo prior distribution parameters and start second time sampling
                try {
                    sameMethylationLevelModel.setModelParameterPseudoPrior();
                    diffMethylationLevelModel.setModelParameterPseudoPrior();
                    this.setModelSamplingList(sameMethylationLevelModel, diffMethylationLevelModel, geneIdx);
                    selectResult = this.sampling(sameMethylationLevelModel, diffMethylationLevelModel, true);

                    // calculate same methylation level model posterior density and diff methylation level model posterior density
                    double diffMethModelProba = this.differentiation(selectResult);
                    double sameMethModelProba = 1 - diffMethModelProba;
                    double bayesFactor = Math.min(diffMethModelProba / sameMethModelProba, 5);
                    // calculate Bayes factor and use the sampling result quantify the methylation level, overdispersion and background expression
                    this.bayesFactors.put(geneIdx, bayesFactor);
                    ModelSelection finalModel = (diffMethModelProba-0.5<0.00001)? sameMethylationLevelModel: diffMethylationLevelModel;
                    double[] finalModelQuantify = this.quantify(finalModel);
                    this.quantifyResult.put(geneIdx, finalModelQuantify);
                    this.modelProbabilities.put(geneIdx, new double[] {sameMethModelProba, diffMethModelProba});

                    // help GC
                    sameMethylationLevelModel.helpGC();
                    diffMethylationLevelModel.helpGC();
                    tretIPReads = null;
                    tretINPUTReads = null;
                    ctrlIPReads = null;
                    ctrlINPUTReads = null;
                    tretIPExpectation = null;
                    tretINPUTExpectation = null;
                    ctrlIPExpectation = null;
                    ctrlINPUTExpectation = null;
                    sameMethylationLevelModel = null;
                    diffMethylationLevelModel = null;
                } catch (Exception e) {
                    e.printStackTrace();
                } finally {
                    countDown.countDown();
                }
            };
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
            sameMethFullConditionProba = sameMethylationLevelModel.iterate(t, sameMethFullConditionProba,modelType==0, usePseudoPrior);
            diffMethFullConditionProba = diffMethylationLevelModel.iterate(t, diffMethFullConditionProba,modelType==1, usePseudoPrior);

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
        // treatment methylation level quantification
        double[] treatmentMethRemainValue = Arrays.stream(model.treatmentMethylationLevel).skip(this.burnIn).sorted().toArray();
        double[] controlMethRemainValue = Arrays.stream(model.controlMethylationLevel).skip(this.burnIn).sorted().toArray();
        double[] treatmentIPOverdispersionRemainValue = Arrays.stream(model.treatmentIPOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] treatmentINPUTOverdispersionRemainValue = Arrays.stream(model.treatmentINPUTOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] controlIPOverdispersionRemainValue = Arrays.stream(model.controlIPOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] controlINPUTOverdispersionRemainValue = Arrays.stream(model.controlINPUTOverdispersion).skip(this.burnIn).sorted().toArray();
        double[] treatmentBkgExpRemainValue = Arrays.stream(model.treatmentBkgExp).skip(this.burnIn).sorted().toArray();
        double[] controlBkgExpRemainValue = Arrays.stream(model.controlBkgExp).skip(this.burnIn).sorted().toArray();

        int medianIdx = remainTime / 2;
        double tretMeth, ctrlMeth, tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion,
               ctrlINPUTOverdispersion, tretBkgExp, ctrlBkgExp;
        if (remainTime%2==0) {
            tretMeth = (treatmentMethRemainValue[medianIdx] + treatmentMethRemainValue[medianIdx+1]) / 2;
            ctrlMeth = (controlMethRemainValue[medianIdx] + controlMethRemainValue[medianIdx+1]) / 2;
            tretIPOverdispersion = (treatmentIPOverdispersionRemainValue[medianIdx] + treatmentIPOverdispersionRemainValue[medianIdx+1]) / 2;
            tretINPUTOverdispersion = (treatmentINPUTOverdispersionRemainValue[medianIdx] + treatmentINPUTOverdispersionRemainValue[medianIdx+1]) / 2;
            ctrlIPOverdispersion = (controlIPOverdispersionRemainValue[medianIdx] + controlIPOverdispersionRemainValue[medianIdx+1]) / 2;
            ctrlINPUTOverdispersion = (controlINPUTOverdispersionRemainValue[medianIdx] + controlINPUTOverdispersionRemainValue[medianIdx+1]) / 2;
            tretBkgExp = (treatmentBkgExpRemainValue[medianIdx] + treatmentBkgExpRemainValue[medianIdx+1]) / 2;
            ctrlBkgExp = (controlBkgExpRemainValue[medianIdx] + controlBkgExpRemainValue[medianIdx+1]) / 2;
        } else {
            tretMeth = treatmentMethRemainValue[medianIdx];
            ctrlMeth = controlMethRemainValue[medianIdx];
            tretIPOverdispersion = treatmentIPOverdispersionRemainValue[medianIdx];
            tretINPUTOverdispersion = treatmentINPUTOverdispersionRemainValue[medianIdx];
            ctrlIPOverdispersion = controlIPOverdispersionRemainValue[medianIdx];
            ctrlINPUTOverdispersion = controlINPUTOverdispersionRemainValue[medianIdx];
            tretBkgExp = treatmentBkgExpRemainValue[medianIdx];
            ctrlBkgExp = controlBkgExpRemainValue[medianIdx];
        }

        return new double[] {tretMeth, ctrlMeth, tretIPOverdispersion, tretINPUTOverdispersion,
                             ctrlIPOverdispersion, ctrlINPUTOverdispersion, tretBkgExp, ctrlBkgExp};
    }

    /**
     * probability of variant methylation level between treatment and control
     * @return probability
     */
    private double differentiation(int[] modelSelectionResult) {
        int remainTime = this.samplingTime - this.burnIn;
        return  (double) (Arrays.stream(modelSelectionResult).skip(this.burnIn).sum()) / remainTime;
    }

    /**
     * treatment and control group IP and INPUT reads count
     * @param geneIdx gene index
     * @return reads count, shape 1 × individualNumber
     */
    private int[][] getGeneReads(int geneIdx) {
        int[] tretIPReads = new int[this.individualNumber], tretINPUTReads = new int[this.individualNumber],
              ctrlIPReads = new int[this.individualNumber], ctrlINPUTReads = new int[this.individualNumber];

        for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
            tretIPReads[individualIdx] = this.treatmentIPReads[individualIdx][geneIdx];
            tretINPUTReads[individualIdx] = this.treatmentINPUTReads[individualIdx][geneIdx];
            ctrlIPReads[individualIdx] = this.controlIPReads[individualIdx][geneIdx];
            ctrlINPUTReads[individualIdx] = this.controlINPUTReads[individualIdx][geneIdx];
        }

        return new int[][] {tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads};
    }

    /**
     * treatment and control group IP and INPUT reads count expectations
     * @param geneIdx gene index
     * @return reads count expectations, shape 1 × individualNumber
     */
    private double[][] getGeneExpectation(int geneIdx) {
        double[] tretIPReads = new double[this.individualNumber], tretINPUTReads = new double[this.individualNumber],
                 ctrlIPReads = new double[this.individualNumber], ctrlINPUTReads = new double[this.individualNumber];

        for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
            tretIPReads[individualIdx] = this.tretIPReadsExpectation[individualIdx][geneIdx];
            tretINPUTReads[individualIdx] = this.tretINPUTReadsExpectation[individualIdx][geneIdx];
            ctrlIPReads[individualIdx] = this.ctrlIPReadsExpectation[individualIdx][geneIdx];
            ctrlINPUTReads[individualIdx] = this.ctrlINPUTReadsExpectation[individualIdx][geneIdx];
        }

        return new double[][] {tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads};
    }

    /**
     * set sampling list for each model
     * @param geneIdx gene index
     */
    private void setModelSamplingList(SameMethylationLevelModel sameMethylationLevelModel, DiffMethylationLevelModel diffMethylationLevelModel, int geneIdx) {
        // here m1, m2 represent same methylation level model and variant methylation level model, respectively
        double[] m1TretMethLevel, m2TretMethLevel, m1TretBkgExp, m2TretBkgExp, m1TretIPOverdispersion, m1TretINPUTOverdispersion, m2TretIPOverdispersion, m2TretINPUTOverdispersion,
                m1CtrlMethLevel, m2CtrlMethLevel, m1CtrlBkgExp, m2CtrlBkgExp, m1CtrlIPOverdispersion, m1CtrlINPUTOverdispersion, m2CtrlIPOverdispersion, m2CtrlINPUTOverdispersion;
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
        sameMethylationLevelModel.setSamplingList(m1TretIPOverdispersion, m1TretINPUTOverdispersion, m1CtrlIPOverdispersion, m1CtrlINPUTOverdispersion, m1TretBkgExp, m1CtrlBkgExp, m1TretMethLevel, m1CtrlMethLevel);
        diffMethylationLevelModel.setSamplingList(m2TretIPOverdispersion, m2TretINPUTOverdispersion, m2CtrlIPOverdispersion, m2CtrlINPUTOverdispersion, m2TretBkgExp, m2CtrlBkgExp, m2TretMethLevel, m2CtrlMethLevel);
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
            bfw.write("#geneIdx\tsameMethModelProba\tdiffMethModelProba\tBF\ttretMeth\tctrlMeth\ttretIPOverdispersion\ttretINPUTOverdispersion\tctrlIPOverdispersion\tctrlINPUTOverdispersion\ttretBkgExp\tctrlBkgExp");
            bfw.newLine();
            String line, records, probabilities;
            double bf;
            double[] quantify, probas;
            for (Integer geneIdx: sortedList) {
                probas = this.modelProbabilities.get(geneIdx);
                probabilities = Arrays.stream(probas).mapToObj(x -> this.df.format((Double) x)).collect(Collectors.joining("\t"));
                quantify = this.quantifyResult.get(geneIdx);
                records = Arrays.stream(quantify).mapToObj(x -> (this.df.format((Double) x))).collect(Collectors.joining("\t"));
                bf = this.bayesFactors.get(geneIdx);
                line = geneIdx + "\t" + probabilities + "\t" + bf + "\t" + records;
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
