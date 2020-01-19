package Quantification;

@FunctionalInterface
public interface SamplingTask {
    Runnable getTask(double prevValues, double curValues, int peakIdx);
}
