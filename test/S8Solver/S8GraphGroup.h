/**
 * S8GraphGroup
 * 
 *
 * Created by Samvel Khalatian on Nov 23, 2010
 * Copyright 2010, All rights reserved
 */

#ifndef S8_GRAPH_GROUP
#define S8_GRAPH_GROUP

#include <map>
#include <memory>
#include <stack>
#include <vector>

#include "S8NumericInput.h"
#include "S8Solution.h"

class TGraph;
class TGraphErrors;
class TObject;

typedef std::vector<NumericInputGroup> BinnedNumericInputGroup;
typedef std::vector<SolutionInBin> BinnedSolution;

struct FlavouredEffGraphGroup
{
    FlavouredEffGraphGroup(const BinnedSolution &);
    FlavouredEffGraphGroup(const BinnedNumericInputGroup &);

    void init(const int &);

    std::auto_ptr<TGraphErrors> b;
    std::auto_ptr<TGraphErrors> cl;
};

struct EffGraphGroup
{
    EffGraphGroup(const std::string &, const BinnedSolution &);
    EffGraphGroup(const std::string &, const BinnedNumericInputGroup &);
    ~EffGraphGroup();

    void draw();
    void save(TDirectory *);

    FlavouredEffGraphGroup mu;
    FlavouredEffGraphGroup tag;

    private:
        typedef std::stack<TObject *> Heap;

        Heap        _heaps;
        std::string _prefix;
};

struct InputGraph
{
    InputGraph(const int &size);

    std::auto_ptr<TGraphErrors> all;
    std::auto_ptr<TGraphErrors> mu;
    std::auto_ptr<TGraphErrors> tag;
    std::auto_ptr<TGraphErrors> muTag;
};

struct InputGraphGroup
{
    InputGraphGroup(const BinnedNumericInputGroup &);
    ~InputGraphGroup();

    void draw();
    void save(TDirectory *);

    InputGraph n;
    InputGraph p;

    private:
        typedef std::stack<TObject *> Heap;

        Heap        _heaps;
};

struct GraphGroup
{
    GraphGroup(const BinnedNumericInputGroup &, const BinnedSolution &);
    ~GraphGroup();

    void save(TDirectory *);
    void draw();

    // Note: objects will be automatically destroyed. Clone graph errors
    // if TMultiGraph is used
    //
    std::auto_ptr<TGraphErrors> alpha;
    std::auto_ptr<TGraphErrors> beta;
    std::auto_ptr<TGraphErrors> gamma;
    std::auto_ptr<TGraphErrors> delta;
    std::auto_ptr<TGraphErrors> kappaB;
    std::auto_ptr<TGraphErrors> kappaCL;

    EffGraphGroup mcEfficiency;
    EffGraphGroup s8Efficiency;

    InputGraphGroup input;

    private:
        typedef std::stack<TObject *> Heap;
        Heap _heaps;
};

#endif