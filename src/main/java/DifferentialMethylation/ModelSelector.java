package DifferentialMethylation;

import org.apache.commons.math3.distribution.UniformRealDistribution;

/**
 * Model selection between two models, the choice based on Bayes factor
 */
public class ModelSelector {
    private UniformRealDistribution urd;

    public ModelSelector() {
        this.urd = new UniformRealDistribution(0, 1);
    }

    /**
     * model selection
     * @return 1, if different methylation level model; otherwise, 0
     */
    public int getSamplingRes(double logSameMethLevelModelPosteriorProba, double logDiffMethLevelModelPosteriorProba) {
//        System.out.println("model1 log proba: " + logSameMethLevelModelPosteriorProba + "\t model2 log proba: " + logDiffMethLevelModelPosteriorProba);
        int modelType = 0;
        double model1Probability = Math.pow(1 + Math.exp(logDiffMethLevelModelPosteriorProba - logSameMethLevelModelPosteriorProba), -1);
//        System.out.println(model1Probability);
        double randNum = this.urd.sample();
        if (randNum > model1Probability)
            return 1;
        return modelType;
    }
}
