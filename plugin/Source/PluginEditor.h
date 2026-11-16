/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class TS808ADAAAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    TS808ADAAAudioProcessorEditor (TS808ADAAAudioProcessor&);
    ~TS808ADAAAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    TS808ADAAAudioProcessor& audioProcessor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TS808ADAAAudioProcessorEditor)
};
