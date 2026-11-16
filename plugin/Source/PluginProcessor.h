/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/**
*/

struct ChannelState
{
    //K-Method
    double w[2] = { 0.0, 0.0 };         // w[n-1]: tensiones de C1 y C2
        double u1 = 0.0;                // u[n-1]
        double y1 = 0.0;                // y[n-1]

    //ADAA
    double u2 = 0.0;                    // u[n-2]: entrada de dos muestras atras
    double pPrev = 0.0;                 // p[n-1]: p de la muestra anterior
    double phPrev[2] = { 0.0, 0.0 };    // ph[n-1]: phat de la muestra anterior

    double FPrev = 0.0;                 // F(p[n-1])
};

struct ModelSet
{
    double Ht[2][2];
	double Gt[2];
    double Jt[2];
    double Kt;
	std::vector<double> pTab, YTab;     // p = x - Kt*y  y su primitiva
};

class TS808ADAAAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    TS808ADAAAudioProcessor();
    ~TS808ADAAAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    //==============================================================================
    juce::AudioProcessorValueTreeState apvts;
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

private:
    //==============================================================================

	void computeModel(double fsEff, ModelSet& m);

    double processK(ChannelState& s, const ModelSet& m, double u);
    double processADAA(ChannelState& s, double u);

	//==============================================================================
	static constexpr double R1 = 4700.0;            // R6
	static constexpr double C1 = 47.0e-9;           // C5
	static constexpr double Rgain = 200000.0;       // R7
	static constexpr double R2 = 51000.0 + Rgain;   // R5 + R7    
	static constexpr double C2 = 51.0e-12;          // C4
	static constexpr double Is = 2.52e-9;           // diodo 1N4148: corriente de saturación
	static constexpr double N = 1.752;              // diodo 1N4148: factor de idealidad
	static constexpr double UT = 25.85e-3;          // diodo 1N4148: voltaje térmico (aprox 27º)

	double A[2][2];                                 // A[fila][columna]
	double B[2];
	double C[2];
	double D[2];

	static constexpr int numPoints = 10001;

	std::vector<double> xTab, yTab;
    
    //==============================================================================
	// prepareToPlay variables

	ModelSet modelK;        // K-method:  fs
	ModelSet modelK2;       // K-method:  2*fs
	ModelSet modelK4;       // K-method:  4*fs
	ModelSet modelADAA;     // ADAA1:     2/3 fs

    ChannelState stateK[2];
    ChannelState stateK2[2];
    ChannelState stateK4[2];
    ChannelState stateADAA[2];

	// LUT variables

    double lookup(const std::vector<double>& pT, const std::vector<double>& vT, double p) const;

    // Oversampling
    
    juce::dsp::Oversampling<float> os2{ 2, 1, juce::dsp::Oversampling<float>::filterHalfBandPolyphaseIIR, true };
    juce::dsp::Oversampling<float> os4{ 2, 2, juce::dsp::Oversampling<float>::filterHalfBandPolyphaseIIR, true };

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TS808ADAAAudioProcessor)   

};
