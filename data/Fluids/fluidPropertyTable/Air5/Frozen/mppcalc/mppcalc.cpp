#include "mutation++.h"
using namespace Mutation;

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

/**
 * Reads a matrix in from a file.
 */
Eigen::MatrixXd readMatrix(std::ifstream& infile)
{
    int cols = 0, rows = 0;
    double buff[((int) 1e6)];

    // Read numbers from file into buffer.
    while (!infile.eof()) {
        std::string line;
        getline(infile, line);
        if (line.empty()) break;

        int temp_cols = 0;
        std::stringstream stream(line);
        while(!stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[cols*i+j];

    return result;
};

/**
 * Base class for all mixture property computers.
 */
class MixtureProperty
{
public:
    typedef Mixture& ARGS;
    static std::string typeName() { return "MixtureProperty"; }
    
    MixtureProperty(ARGS mix, string name, string doc, string units) : 
        m_mix_ref(mix),
        m_name(name),
        m_doc_string(doc),
        m_units(units)
    { } 
    
    virtual ~MixtureProperty() {}
    
    // Getters
    const string& name() const { return m_name; }
    const string& docString() const { return m_doc_string; }
    const string& units() const { return m_units; }
    
    virtual void compute(VectorXd& v, VectorXd& args) { };
    //virtual void compute(VectorXd& v) { }

protected:

    Mixture& mix() const { return m_mix_ref; }

private:
    
    Mixture& m_mix_ref;
    std::string m_name;
    std::string m_doc_string;
    std::string m_units;
};


/**
 * Property parser
 */
class PropertyParser
{
public:
    PropertyParser(const std::string& to_parse, Mixture& mix)
    {
        // Need to extract the property name and the arguement indices
        std::string name = to_parse;
        std::string indices;

        // Look for (...) at end of name
        std::size_t lpi = name.find_last_of('(');
        if (lpi != std::string::npos && name[name.size()-1] == ')') {
            name = to_parse.substr(0, lpi);
            indices = to_parse.substr(lpi+1, to_parse.size()-lpi-2);
        }

        // Load the property
        m_property = shared_ptr<MixtureProperty>(
            Utilities::Config::Factory<MixtureProperty>::create(name, mix));

        // Now load the indices
        if (!parseIndices(indices)) {
            cout << "Incorrect indices for property " << name << "!" << endl;
            exit(1);
        }

        m_args.resize(m_args_indices.size());
    }
    
    shared_ptr<MixtureProperty> property() const { return m_property; }

    template <typename T>
    const VectorXd& compute(const DenseBase<T>& row) {
        // Copy over arguments (if there are any)
        if (m_args_indices.size() > 0)
            for (int i = 0; i < m_args_indices.size(); ++i)
                m_args[i] = row[m_args_indices[i]];

        // Compute the property and return the result
        m_property->compute(m_result, m_args);
        return m_result;
    };

private:

    bool parseIndices(const std::string& list)
    {
        std::vector<std::string> ranges;
        std::vector<std::string> bounds;
        Utilities::String::tokenize(list, ranges, ",");

        std::vector<std::string>::const_iterator iter = ranges.begin();
        for ( ; iter != ranges.end(); ++iter) {
            bounds.clear();
            Utilities::String::tokenize(*iter, bounds, "-");

            if (!Utilities::String::isNumeric(bounds))
                return false;

            switch (bounds.size()) {
                case 1: {
                    int i = atoi(bounds[0].c_str());
                    if (i < 0)
                        return false;
                    m_args_indices.push_back(atoi(bounds[0].c_str())-1);
                    break;
                }
                case 2: {
                    int i1 = atoi(bounds[0].c_str()-1);
                    int i2 = atoi(bounds[1].c_str()-1);

                    if (i1 >= i2 || i1 < 0)
                        return false;

                    for (int i = i1; i <= i2; ++i)
                        m_args_indices.push_back(i);

                    break;
                }
                default:
                    return false;
            }
        }

        return true;
    }

private:
    shared_ptr<MixtureProperty> m_property;
    std::vector<int> m_args_indices;

    VectorXd m_args;
    VectorXd m_result;
};


// Macro for adding a new mixture property
#define ADD_MIXTURE_PROPERTY(__NAME__,__DOCS__,__UNITS__,__SIZE__,__CODE__)\
class __NAME__ : public MixtureProperty\
{\
public:\
    __NAME__(Mixture& mix) :\
        MixtureProperty(mix, #__NAME__, __DOCS__, __UNITS__) {}\
    void compute(VectorXd& v, VectorXd& args) {\
        v.resize(__SIZE__);\
        __CODE__;\
    };\
};\
Utilities::Config::ObjectProvider<\
    __NAME__, MixtureProperty> op_##__NAME__(#__NAME__);
    

// Densities
ADD_MIXTURE_PROPERTY(rho,
    "mixture density", "kg/m^3", 1, 
    v[0] = mix().density()
)

ADD_MIXTURE_PROPERTY(rhoi,
    "species densities", "kg/m^3", mix().nSpecies(), 
    v = Map<const VectorXd>(mix().Y(), mix().nSpecies()) * mix().density()
)

// Mass Fractions
ADD_MIXTURE_PROPERTY(yi,
    "species mass fractions", "", mix().nSpecies(), 
    v = Map<const VectorXd>(mix().Y(), mix().nSpecies())
)

// Mole Fractions
ADD_MIXTURE_PROPERTY(xi,
    "species mole fractions", "", mix().nSpecies(), 
    v = Map<const VectorXd>(mix().X(), mix().nSpecies())
)

// Pressures
ADD_MIXTURE_PROPERTY(p,
    "mixture pressure", "Pa/m^2", 1, 
    v[0] = mix().P()
)

// Production rates
ADD_MIXTURE_PROPERTY(omegai,
    "species production rates due to chemical reactions", "kg/m^3-s", mix().nSpecies(),
    mix().netProductionRates(v.data())
)

// Thermodynamics
ADD_MIXTURE_PROPERTY(cp, 
    "mixture specific heat at constant pressure", "J/kg-K", 1,
    v[0] = mix().mixtureFrozenCpMass()
)

ADD_MIXTURE_PROPERTY(cpi, 
    "species specific heat at constant pressure", "J/kg-K", mix().nSpecies(),
    mix().speciesCpOverR(v.data());
    for (int i = 0; i < mix().nSpecies(); ++i)
        v[i] *= (RU / mix().speciesMw(i));
)

ADD_MIXTURE_PROPERTY(Cpi, 
    "species specific heat at constant pressure", "J/mol-K", mix().nSpecies(),
    mix().speciesCpOverR(v.data());
    v *= RU;
)

ADD_MIXTURE_PROPERTY(h, 
    "mixture enthalpy", "J/kg", 1,
    v[0] = mix().mixtureHMass()
)

ADD_MIXTURE_PROPERTY(hi, 
    "species enthalpies", "J/kg", mix().nSpecies(),
    mix().speciesHOverRT(v.data());
    for (int i = 0; i < mix().nSpecies(); ++i)
        v[i] *= (RU * mix().T() / mix().speciesMw(i));
)

ADD_MIXTURE_PROPERTY(H, 
    "mixture enthalpy", "J/mol", 1,
    v[0] = mix().mixtureHMole()
)

ADD_MIXTURE_PROPERTY(Hi, 
    "species enthalpies", "J/mol", mix().nSpecies(),
    mix().speciesHOverRT(v.data());
    v *= (RU * mix().T());
)

ADD_MIXTURE_PROPERTY(si, 
    "species entropies", "J/kg-K", mix().nSpecies(),
    mix().speciesSOverR(v.data());
    for (int i = 0; i < mix().nSpecies(); ++i)
        v[i] *= (RU / mix().speciesMw(i));
)

ADD_MIXTURE_PROPERTY(Si, 
    "species entropies", "J/mol-K", mix().nSpecies(),
    mix().speciesSOverR(v.data());
    v *= RU;
)

// Temperatures
ADD_MIXTURE_PROPERTY(Th,
    "heavy particle translational temperature", "K", 1,
    v(0) = mix().T())

ADD_MIXTURE_PROPERTY(Te,
    "free electron translational temperature", "K", 1,
    v(0) = mix().Te())
    
// Thermal conductivities
ADD_MIXTURE_PROPERTY(lambda, 
    "frozen thermal conductivity", "W/m-K", 1,
    v(0) = mix().frozenThermalConductivity())

ADD_MIXTURE_PROPERTY(lambda_eq, 
    "equilibrium thermal conductivity", "W/m-K", 1,
    v(0) = mix().equilibriumThermalConductivity())

ADD_MIXTURE_PROPERTY(lambda_h, 
    "heavy particle thermal conductivity", "W/m-K", 1,
    v(0) = mix().heavyThermalConductivity())

ADD_MIXTURE_PROPERTY(lambda_e, 
    "electron thermal conductivity", "W/m-K", 1,
    v(0) = mix().electronThermalConductivity())
    
ADD_MIXTURE_PROPERTY(lambda_rot, 
    "rotational energy thermal conductivity", "W/m-K", 1,
    v(0) = mix().rotationalThermalConductivity())

ADD_MIXTURE_PROPERTY(lambda_vib, 
    "internal energy thermal conductivity", "W/m-K", 1,
    v(0) = mix().vibrationalThermalConductivity())

ADD_MIXTURE_PROPERTY(lambda_elec, 
    "electronic energy thermal conductivity", "W/m-K", 1,
    v(0) = mix().electronicThermalConductivity())

// Viscosity
ADD_MIXTURE_PROPERTY(mu, 
    "dynamic viscosity", "Pa-s", 1,
    v(0) = mix().viscosity())

// Fake "properties"
ADD_MIXTURE_PROPERTY(col,
    "prints columns in parentheses", "", args.size(),
    v = args)



void printAvailableProps()
{
    // First get a list of all the property names that are registered
    vector<string> names = Utilities::Config::Factory<MixtureProperty>::names();
    
    // Make sure the list is ordered
    std::sort(names.begin(), names.end());
    
    // Loop over each name and return the required data
    Mixture mix("air5");
    MixtureProperty* p_prop;
    for (int i = 0; i < names.size(); ++i) {
        p_prop = Utilities::Config::Factory<MixtureProperty>::create(
            names[i], mix);
        cout << setw(20) << left << p_prop->name() 
             << setw(15) << left << p_prop->units()
             << p_prop->docString()
             << endl;
        delete p_prop;
    }
}

/**
 * Driver program for comparing a function to the values stored in the given
 * file.
 */
int main(int argc, char* argv[])
{
    if (argc == 1) {
        cout << "\nAvailable properties: " << endl;
        printAvailableProps();
        cout << endl;
        return 0;
    }
    
    // Open the file
    ifstream file(argv[1]);

    // Load the mixture
    string str;
    getline(file, str);
    MixtureOptions opts(str);
    getline(file, str);
    opts.setStateModel(str);
    Mixture mix(opts);

    // Load the state variable information
    getline(file, str);
    int var_set, s1, n1, s2, n2;
    stringstream ss(str);  ss >> var_set >> s1 >> n1 >> s2 >> n2;

    // Load the functions to compute
    getline(file, str);

    std::vector<std::string> function_names;
    Utilities::String::tokenize(str, function_names, " ");
    
    std::vector<PropertyParser> props;
    for (int i = 0; i < function_names.size(); ++i)
        props.push_back(PropertyParser(function_names[i], mix));

    // Load the state variables and data to compare to
    MatrixXd data = readMatrix(file);

    // Close the file
    file.close();

    // Loop over the states
    double v1 [n1];
    double v2 [n2];
    double args [data.cols()];
    
    bool matches = true;

    for (int i = 0; i < data.rows(); ++i) {
        // Set the state of the mixture
        Map<VectorXd>(v1,n1) = data.row(i).segment(s1-1,n1);
        Map<VectorXd>(v2,n2) = data.row(i).segment(s2-1,n2);
        mix.setState(v1, v2, var_set);

        // Now compute each function
        for (int j = 0; j < props.size(); ++j) {
            const VectorXd& res = props[j].compute(data.row(i));
            for (int k = 0; k < res.size(); ++k)
                cout << setw(15) << res[k] << " ";
        }
        cout << endl;
    }

    return 0;
}
