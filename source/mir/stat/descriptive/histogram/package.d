/++
This module contains algorithms for create histograms.

License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

Authors: John Michael Hall, Ilya Yaroshenko

Copyright: 2020 Mir Stat Authors.

+/

module mir.stat.descriptive.histogram;

public import mir.stat.descriptive.histogram.frequency;
public import mir.stat.descriptive.histogram.breaks;
public import mir.stat.descriptive.histogram.accumulator;
public import mir.stat.descriptive.histogram.axis;
public import mir.stat.descriptive.histogram.traits;
public import mir.stat.descriptive.histogram.api;

//date: 5/17/2021
// helper functions in axis need documented UTs
// got transformAxis working properly, need to add documentation for TODOs, 
// rchistogram is missing UTs for transformAxis new version

//later TODOs
//see if it is possible to get multiple axis working?
//use binarySearch for search instead of loop for variable
//provide frequency
//provide cumulative frequency function

//frequency can be implemented separately from histogram by handling counts in that

// 11/12/2020
// Priorities: 
// 1) a) Add slice functions
//    b) Add range interface
// 2) A way to convert the results to string that Ilya will be happy with
// 3) Multi-axis histograms
// 4) Historgram with GC version, makeHistogram to handle any allocation strategy

// Not part of MVP
// 1) Make Storage able handle other ways to count, such as  MeanAccumulator, like dense_storage for boost histogram
// 2) FrequencyAccumulator/Frequency/CumFrequency
// 3) weighted histogram/frequency
// 4) Replace phobos in VariableAxis
// 5) Better count type (can increase itself)

