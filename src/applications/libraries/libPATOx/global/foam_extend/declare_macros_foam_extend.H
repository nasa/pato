// Foam Extend 4.1: macros

// PATO
#define lookup_foam_extend_mesh(mesh_foam,region_name)   \
const_cast<Foam_extend_::fvMesh&>(                       \
mesh_foam.time().db().lookupObject<simpleMaterialsModel> \
("MaterialsModel").get_foam_extend_mesh(region_name))


// Foam_extend_fvPatchField.H
#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)   \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, patch                                 \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField,                                                           \
    typePatchTypeField,                                                       \
    patchMapper                                                               \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, dictionary                            \
);


#define makePatchTypeFieldTypeName(type)                                      \
                                                                              \
defineNamedTemplateTypeNameAndDebug(type, 0);


#define makePatchFieldsTypeName(type)                                         \
                                                                              \
makePatchTypeFieldTypeName(type##FvPatchScalarField);                         \
makePatchTypeFieldTypeName(type##FvPatchVectorField);                         \
makePatchTypeFieldTypeName(type##FvPatchSphericalTensorField);                \
makePatchTypeFieldTypeName(type##FvPatchSymmTensorField);                     \
makePatchTypeFieldTypeName(type##FvPatchSymmTensor4thOrderField);             \
makePatchTypeFieldTypeName(type##FvPatchDiagTensorField);                     \
makePatchTypeFieldTypeName(type##FvPatchTensorField);


#define makePatchTypeField(PatchTypeField, typePatchTypeField)                \
                                                                              \
defineTypeNameAndDebug(typePatchTypeField, 0);                                \
                                                                              \
addToPatchFieldRunTimeSelection                                               \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);


#define makeTemplatePatchTypeField(PatchTypeField, typePatchTypeField)        \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                   \
                                                                              \
addToPatchFieldRunTimeSelection                                               \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);


#define makePatchFields(type)                                                 \
                                                                              \
makeTemplatePatchTypeField(fvPatchScalarField, type##FvPatchScalarField);     \
makeTemplatePatchTypeField(fvPatchVectorField, type##FvPatchVectorField);     \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSphericalTensorField,                                              \
    type##FvPatchSphericalTensorField                                         \
);                                                                            \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSymmTensorField,                                                   \
    type##FvPatchSymmTensorField                                              \
);                                                                            \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSymmTensor4thOrderField,                                           \
    type##FvPatchSymmTensor4thOrderField                                      \
);                                                                            \
makeTemplatePatchTypeField                                                    \
 (                                                                            \
    fvPatchDiagTensorField,                                                   \
    type##FvPatchDiagTensorField                                              \
);                                                                            \
makeTemplatePatchTypeField(fvPatchTensorField, type##FvPatchTensorField);


#define makePatchTypeFieldTypedefs(type)                                      \
                                                                              \
typedef type##FvPatchField<scalar> type##FvPatchScalarField;                  \
typedef type##FvPatchField<vector> type##FvPatchVectorField;                  \
typedef type##FvPatchField<sphericalTensor>                                   \
    type##FvPatchSphericalTensorField;                                        \
typedef type##FvPatchField<symmTensor> type##FvPatchSymmTensorField;          \
typedef type##FvPatchField<symmTensor4thOrder>                                \
    type##FvPatchSymmTensor4thOrderField;                                     \
typedef type##FvPatchField<diagTensor>                                        \
    type##FvPatchDiagTensorField;                                             \
typedef type##FvPatchField<tensor> type##FvPatchTensorField;


// Foam_extend_error.H
//- Report an error message using Foam_extend_::FatalError
//  for functionName in file __FILE__ at line __LINE__
#define FatalErrorIn(functionName)                                             \
    ::Foam_extend_::FatalError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam_extend_::FatalError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)


//- Report an error message using Foam_extend_::FatalIOError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorIn(functionName, ios)                                      \
    ::Foam_extend_::FatalIOError((functionName), __FILE__, __LINE__, (ios))

//- Report an error message using Foam_extend_::FatalIOError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorInFunction(ios) FatalIOErrorIn(FUNCTION_NAME, ios)


//- Report an error message using Foam_extend_::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorIn(functionName, ios, msg)                             \
    ::Foam_extend_::IOerror::SafeFatalIOError                                          \
    ((functionName), __FILE__, __LINE__, (ios), (msg))

//- Report an error message using Foam_extend_::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorInFunction(ios, msg)                                   \
    SafeFatalIOErrorIn(FUNCTION_NAME, ios, msg)


//- Issue a FatalErrorIn for a function not currently implemented.
//  The functionName is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define notImplemented(functionName)                                           \
    FatalErrorIn(functionName)                                                 \
        << "Not implemented" << ::Foam_extend_::abort(FatalError);

//- Issue a FatalErrorIn for a function not currently implemented.
//  The FUNCTION_NAME is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define NotImplemented notImplemented(FUNCTION_NAME)


// No defineDebugSwitch.H

// Foam_extend_className.H
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declarations (without debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//- Add typeName information from argument @a TypeNameString to a class.
//  Without debug information
#define ClassNameNoDebug(TypeNameString)                                      \
    static const char* typeName_() { return TypeNameString; }                 \
    static const ::Foam_extend_::word typeName

//- Add typeName information from argument @a TypeNameString to a namespace.
//  Without debug information.
#define NamespaceNameNoDebug(TypeNameString)                                  \
    inline const char* typeName_() { return TypeNameString; }                 \
    extern const ::Foam_extend_::word typeName

//- Add typeName information from argument @a TemplateNameString to a template class.
//  Without debug information.
#define TemplateNameNoDebug(TemplateNameString)                               \
class TemplateNameString##Name                                                \
{                                                                             \
public:                                                                       \
    TemplateNameString##Name() {}                                             \
    ClassNameNoDebug(#TemplateNameString);                                    \
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declarations (with debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Add typeName information from argument @a TypeNameString to a class.
//  Also declares debug information.
#define ClassName(TypeNameString)                                             \
    ClassNameNoDebug(TypeNameString);                                         \
    static Foam_extend_::debug::debugSwitch debug;

//- Add typeName information from argument @a TypeNameString to a namespace.
//  Also declares debug information.
#define NamespaceName(TypeNameString)                                         \
    NamespaceNameNoDebug(TypeNameString);                                     \
    extern Foam_extend_::debug::debugSwitch debug;

//- Add typeName information from argument @a TypeNameString to a template class.
//  Also declares debug information.
#define TemplateName(TemplateNameString)                                      \
class TemplateNameString##Name                                                \
{                                                                             \
public:                                                                       \
    TemplateNameString##Name() {}                                             \
    ClassName(#TemplateNameString);                                           \
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (without debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the typeName, with alternative lookup as @a Name
#define defineTypeNameWithName(Type, Name)                                    \
    const ::Foam_extend_::word Type::typeName(Name)

//- Define the typeName
#define defineTypeName(Type)                                                  \
    defineTypeNameWithName(Type, Type::typeName_())

#ifdef __INTEL_COMPILER
//- Define the typeName as @a Name for template classes
# define defineTemplateTypeNameWithName(Type, Name)                           \
    defineTypeNameWithName(Type, Name)
#else
//- Define the typeName as @a Name for template classes
# define defineTemplateTypeNameWithName(Type, Name)                           \
    template<>                                                                \
    defineTypeNameWithName(Type, Name)
#endif

//- Define the typeName for template classes, useful with typedefs
#define defineTemplateTypeName(Type)                                          \
    defineTemplateTypeNameWithName(Type, #Type)

//- Define the typeName directly for template classes
#define defineNamedTemplateTypeName(Type)                                     \
    defineTemplateTypeNameWithName(Type, Type::typeName_())



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (debug information only)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the debug information, lookup as @a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch, SwitchDescr)       \
    ::Foam_extend_::debug::debugSwitch                                                \
    Type::debug(::Foam_extend_::word(Name), DebugSwitch, SwitchDescr)

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch, SwitchDescr)                     \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch, SwitchDescr);

#ifdef __INTEL_COMPILER
//- Define the debug information for templates, lookup as @a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    defineDebugSwitchWithName(Type, Name, DebugSwitch, "");

#else
//- Define the debug information for templates, lookup as @a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    template<>                                                                \
    defineDebugSwitchWithName(Type, Name, DebugSwitch, "");
#endif

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, DebugSwitch)                          \
    defineTemplateDebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, DebugSwitch)                     \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (with debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the typeName and debug information
#define defineTypeNameAndDebug(Type, DebugSwitch)                             \
    defineTypeName(Type);                                                     \
    defineDebugSwitch(Type, DebugSwitch, "")

//- Define the typeName and debug information + description
#define defineTypeNameAndDebugWithDescription(Type, DebugSwitch, SwitchDescr) \
    defineTypeName(Type);                                                     \
    defineDebugSwitch(Type, DebugSwitch, SwitchDescr)

//- Define the typeName and debug information, lookup as @a Name
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)       \
    defineTemplateTypeNameWithName(Type, Name);                               \
    defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)

//- Define the typeName and debug information for templates, useful with typedefs
#define defineTemplateTypeNameAndDebug(Type, DebugSwitch)                     \
    defineTemplateTypeNameAndDebugWithName(Type, #Type, DebugSwitch)

//- Define the typeName and debug information for templates
#define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                \
    defineNamedTemplateTypeName(Type);                                        \
    defineNamedTemplateDebugSwitch(Type, DebugSwitch)

