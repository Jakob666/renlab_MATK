package DifferentialMethylation;

import Quantification.BackgroundExpressionSampler;

@FunctionalInterface
public interface BkgExpSamplingTask {
    Runnable createSamplingTask(double prevValue, double curValue, BackgroundExpressionSampler sampler, int idx);
}
