/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
TS808ADAAAudioProcessor::TS808ADAAAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
    : AudioProcessor(BusesProperties()
#if ! JucePlugin_IsMidiEffect
#if ! JucePlugin_IsSynth
        .withInput("Input", juce::AudioChannelSet::stereo(), true)
#endif
        .withOutput("Output", juce::AudioChannelSet::stereo(), true)
#endif
    ), apvts(*this, nullptr, "State", createParameterLayout())
#endif
{

    //Matrix A
    A[0][0] = -1.0 / (R1 * C1);
    A[0][1] = 0;
    A[1][0] = -1.0 / (R1 * C2);
    A[1][1] = -1.0 / (R2 * C2);

    //Matrix B
    B[0] = 1.0 / (R1 * C1);
    B[1] = 1.0 / (R1 * C2);

    //Matrix C
    C[0] = 0;
    C[1] = -1.0 / C2;

    //Matrix D
    D[0] = 0;
    D[1] = 1;

    // Curva del diodo y = f(x), de 0 a 1 V

    xTab.resize(numPoints);
    yTab.resize(numPoints);

    // K-method
    modelK.pTab.resize(numPoints);
    modelK.YTab.resize(numPoints);

    // K-method 2x
	modelK2.pTab.resize(numPoints);
	modelK2.YTab.resize(numPoints);

    // K-method 4x
	modelK4.pTab.resize(numPoints);
	modelK4.YTab.resize(numPoints);

	// ADAA1
    modelADAA.pTab.resize(numPoints);
    modelADAA.YTab.resize(numPoints);

    for (int i = 0; i < numPoints; ++i) {

        xTab[i] = i / 10000.0;                                                                  // x values from 0 to 1
        yTab[i] = Is * (std::exp(xTab[i] / (N * UT)) - std::exp(-xTab[i] / (N * UT)));         // y = Is * (exp(x/(N*UT)) - exp(-x/(N*UT)))
    }

    DBG("Constructor | A11 = " << A[1][1] << " | y(0.5 V) = " << yTab[5000]);
}

TS808ADAAAudioProcessor::~TS808ADAAAudioProcessor()
{
}

//==============================================================================
const juce::String TS808ADAAAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool TS808ADAAAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}

bool TS808ADAAAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}

bool TS808ADAAAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
}

double TS808ADAAAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int TS808ADAAAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
    // so this should be at least 1, even if you're not really implementing programs.
}

int TS808ADAAAudioProcessor::getCurrentProgram()
{
    return 0;
}

void TS808ADAAAudioProcessor::setCurrentProgram(int index)
{
}

const juce::String TS808ADAAAudioProcessor::getProgramName(int index)
{
    return {};
}

void TS808ADAAAudioProcessor::changeProgramName(int index, const juce::String& newName)
{
}

void TS808ADAAAudioProcessor::computeModel(double fsEff, ModelSet& m)
{
    const double h = 2.0 * fsEff;

    // M = h*I - A y su inversa
    double a = h - A[0][0];
    double b = -A[0][1];
    double c = -A[1][0];
    double d = h - A[1][1];

    double detM = a * d - b * c;
    double invA[2][2] = { {d / detM, -b / detM}, {-c / detM, a / detM} };

    // Ht = invA * (h*I + A)
    double P[2][2] = { {h + A[0][0], A[0][1]}, {A[1][0], h + A[1][1]} };
    m.Ht[0][0] = invA[0][0] * P[0][0] + invA[0][1] * P[1][0];
    m.Ht[0][1] = invA[0][0] * P[0][1] + invA[0][1] * P[1][1];
    m.Ht[1][0] = invA[1][0] * P[0][0] + invA[1][1] * P[1][0];
    m.Ht[1][1] = invA[1][0] * P[0][1] + invA[1][1] * P[1][1];

    // Gt = invA * C
    m.Gt[0] = invA[0][0] * C[0] + invA[0][1] * C[1];
    m.Gt[1] = invA[1][0] * C[0] + invA[1][1] * C[1];

    // Jt = invA * B
    m.Jt[0] = invA[0][0] * B[0] + invA[0][1] * B[1];
    m.Jt[1] = invA[1][0] * B[0] + invA[1][1] * B[1];

    // Kt = D * Gt
    m.Kt = D[0] * m.Gt[0] + D[1] * m.Gt[1];

    //LUT
    m.YTab[0] = 0;
    m.pTab[0] = xTab[0] - m.Kt * yTab[0];

    for (int i = 1; i < numPoints; ++i) {
        m.pTab[i] = xTab[i] - m.Kt * yTab[i];
        m.YTab[i] = m.YTab[i - 1] + 0.5 * (yTab[i] + yTab[i - 1]) * (m.pTab[i] - m.pTab[i - 1]);
    }
}

double TS808ADAAAudioProcessor::lookup(const std::vector<double>& pT, const std::vector<double>& vT, double p) const
{
    // Aplicamos clamping para valores fuera del rango de la tabla
    if (p <= pT[0]) return vT[0];
    else if (p >= pT[numPoints - 1]) return vT[numPoints - 1];

    // Búsqueda binaria para encontrar el índice correspondiente
    int low = 0;
    int high = numPoints - 1;
    int mid;
    while (high - low > 1) {
        mid = (low + high) / 2;
        if (pT[mid] <= p) low = mid;
        else high = mid;
    }

	// Interpolación lineal entre los dos puntos encontrados
    double t = (p - pT[low]) / (pT[low + 1] - pT[low]);         // Entre 0 (p = p0) y 1 (p = p1)
    
    return vT[low] + (vT[low + 1] - vT[low]) * t;
}

double TS808ADAAAudioProcessor::processK(ChannelState& s, const ModelSet& m, double u)
{
	// pasos 1 a 6 del K-method:

		// Paso 1: p_hat = Ht*w[n-1] + Jt*(u[n] + u[n-1]) + Gt*y[n-1]
		// Todo lo que ya se conoce de esta muestra (la parte lineal), sin la corriente actual del diodo.
		// p_hat tiene 2 componentes: ph0 (fila 0) y ph1 (fila 1).

	const double ph0 = (m.Ht[0][0] * s.w[0] + m.Ht[0][1] * s.w[1])
		+ m.Jt[0] * (u + s.u1) + m.Gt[0] * s.y1;
	const double ph1 = (m.Ht[1][0] * s.w[0] + m.Ht[1][1] * s.w[1])
		+ m.Jt[1] * (u + s.u1) + m.Gt[1] * s.y1;

	// Paso 2: tension que "ve" el diodo
	const double p = ph1;          // p = D * p_hat, con D = [0 1]

	// Paso 3: corriente del diodo (parte no lineal) desde la LUT del K-method.
	// La tabla solo tiene la mitad positiva: se busca |p| y se recupera el signo (diodos simetricos)
	double y = lookup(m.pTab, yTab, std::abs(p));
	if (p < 0) y = -y;

	// Paso 4: nuevas tensiones de los condensadores, w[n] = Gt*y[n] + p_hat
	const double wNew[2] = { m.Gt[0] * y + ph0, m.Gt[1] * y + ph1 };

	// Paso 5: salida de la etapa de recorte = tension de C2 + entrada
	const double out = wNew[1] + u;

	// Paso 6: guardar la memoria para la siguiente muestra (siempre al final)
	s.w[0] = wNew[0];
	s.w[1] = wNew[1];
	s.u1 = u;
	s.y1 = y;

	return out;
}

double TS808ADAAAudioProcessor::processADAA(ChannelState& s, double u) {

    const double ph0 = (modelADAA.Ht[0][0] * s.w[0] + modelADAA.Ht[0][1] * s.w[1])
        + modelADAA.Jt[0] * (u + 0.5 * (s.u1 + s.u2)) + modelADAA.Gt[0] * s.y1;
    const double ph1 = (modelADAA.Ht[1][0] * s.w[0] + modelADAA.Ht[1][1] * s.w[1])
        + modelADAA.Jt[1] * (u + 0.5 * (s.u1 + s.u2)) + modelADAA.Gt[1] * s.y1;

    const double p = ph1;

    const double dp = p - s.pPrev;

    double y = 0.0;

    const double tol = std::numeric_limits<double>::epsilon();   // 2.2e-16, igual que eps en MATLAB
    const double F0 = lookup(modelADAA.pTab, modelADAA.YTab, std::abs(p));
    if (std::abs(dp) < tol) {
        double mid = 0.5 * (p + s.pPrev);
        y = lookup(modelADAA.pTab, yTab, std::abs(mid));
        if (mid < 0) y = -y;
    }
    else {
        y = (F0 - s.FPrev) / (dp);
    }

    const double wNew[2] = { modelADAA.Gt[0] * y + 0.5 * (ph0 + s.phPrev[0]), modelADAA.Gt[1] * y + 0.5 * (ph1 + s.phPrev[1]) };
    s.w[0] = wNew[0];
    s.w[1] = wNew[1];
    const double out = wNew[1] + u;
    s.u2 = s.u1;
    s.u1 = u;
    s.y1 = y;
    s.pPrev = p;
    s.phPrev[0] = ph0;
    s.phPrev[1] = ph1;
	s.FPrev = F0;

    return out;
}

//==============================================================================
void TS808ADAAAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
	computeModel(sampleRate,                modelK);            // K-method: fs
	computeModel(sampleRate * 2.0,          modelK2);           // K-method 2x: fs*2
    computeModel(sampleRate * 4.0,          modelK4);           // K-method 4x: fs*4
    computeModel(sampleRate * 2.0 / 3.0,    modelADAA);         // ADAA-method: fs*2/3
    
    os2.initProcessing(samplesPerBlock);
    os2.reset();

    os4.initProcessing(samplesPerBlock);
    os4.reset();

    for (int ch = 0; ch < 2; ++ch) {
		stateK[ch] = ChannelState();
		stateK2[ch] = ChannelState();
		stateK4[ch] = ChannelState();   
        stateADAA[ch] = ChannelState();
    }

    DBG("prepareToPlay | fs = " << sampleRate << " Hz");
    DBG("modelK    | Kt = " << modelK.Kt << " | p[5000] = " << modelK.pTab[5000] << " | Y[10000] = " << modelK.YTab[10000]);
    DBG("modelADAA | Kt = " << modelADAA.Kt << " | p[5000] = " << modelADAA.pTab[5000] << " | Y[10000] = " << modelADAA.YTab[10000]);
    DBG("lookup | fila 5000 = " << lookup(modelK.pTab, yTab, modelK.pTab[5000]) << " | p=0 -> " << lookup(modelK.pTab, yTab, 0.0) << " | p=1e9 -> " << lookup(modelK.pTab, yTab, 1e9));
}

void TS808ADAAAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool TS808ADAAAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused(layouts);
    return true;
#else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
#endif

    return true;
#endif
}
#endif

void TS808ADAAAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    const int method = juce::roundToInt(apvts.getRawParameterValue("method")->load());


    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear(i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.

    // Parametros: se leen una vez por bloque y se pasan de dB a factor lineal
    const float inputDb = apvts.getRawParameterValue("input")->load();
    const float inGain = juce::Decibels::decibelsToGain(inputDb);

    const float outputDb = apvts.getRawParameterValue("output")->load();
    const float outGain = juce::Decibels::decibelsToGain(outputDb);

	// Un "pedal virtual" por canal: cada canal (L, R) usa su propia memoria s
    if (method == 0 || method == 1) {
        for (int channel = 0; channel < totalNumInputChannels; ++channel) {
            auto& sK = stateK[channel];                         // memoria del K-method de este canal
            auto& sADAA = stateADAA[channel];                   // memoria del ADAA1 de este canal

            auto* channelData = buffer.getWritePointer(channel);

            for (int i = 0; i < buffer.getNumSamples(); i++) {
                const double u = (double)channelData[i] * inGain;
                double out = 0.0;

                switch (method)
                {
                case 0:
                    out = processK(sK, modelK, u);
                    break;
                case 1:
                    out = processADAA(sADAA, u);
                    break;
                default:
                    break;
                }

                channelData[i] = (float)(out * outGain);
            }
        }
    }
    else {
		const bool isx4 = (method == 3);

        
        auto& os = isx4 ? os4 : os2;
        auto& m = isx4 ? modelK4 : modelK2;
		auto* s = isx4 ? stateK4 : stateK2;

        juce::dsp::AudioBlock<float> block(buffer);
        auto up = os.processSamplesUp(block);

        for (int channel = 0; channel < totalNumInputChannels; ++channel) {
            auto* data = up.getChannelPointer(channel);
            for (int i = 0; i < (int)up.getNumSamples(); i++) {
                const double u = (double)data[i] * inGain;
                const double out = processK(s[channel], m, u);

				data[i] = (float)(out * outGain);
            }
        }

        os.processSamplesDown(block);

    }

}

//==============================================================================
bool TS808ADAAAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* TS808ADAAAudioProcessor::createEditor()
{
    return new juce::GenericAudioProcessorEditor(*this);
}

//==============================================================================
void TS808ADAAAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.

    auto state = apvts.copyState();
    auto stateXml = state.createXml();

    copyXmlToBinary(*stateXml, destData);
}

void TS808ADAAAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.

    auto xmlState = getXmlFromBinary(data, sizeInBytes);
    if (xmlState != nullptr && xmlState->hasTagName(apvts.state.getType())) {
        auto tree = juce::ValueTree::fromXml(*xmlState);
        apvts.replaceState(tree);
    }
}

juce::AudioProcessorValueTreeState::ParameterLayout
TS808ADAAAudioProcessor::createParameterLayout() {
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> p;
    p.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID{ "input", 1 }, "Input (dB)", 
        -24.0f /* min */, 24.0f /* max */, 0.0f /* defecto */));

    p.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID{ "output", 1 }, "Output (dB)", 
        -24.0f /* min */, 0.0f /* max */, -6.0f /* defecto */));

    p.push_back(std::make_unique<juce::AudioParameterChoice>(
        juce::ParameterID{ "method", 1 }, "Method", juce::StringArray{ "K-method", "ADAA1", "K x2", "K x4"}, 0));
    
    return{ p.begin(), p.end() };
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TS808ADAAAudioProcessor();
}