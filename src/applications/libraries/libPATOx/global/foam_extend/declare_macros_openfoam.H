// OpenFOAM 7: macros

// fvPatchField.H
#define makeFvPatchField(fvPatchTypeField)                                     \
                                                                               \
defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                      \
template<>                                                                     \
int fvPatchTypeField::disallowGenericFvPatchField                              \
(                                                                              \
    debug::debugSwitch("disallowGenericFvPatchField", 0)                       \
);                                                                             \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);                  \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);            \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);


#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)    \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// Use with caution
#define addRemovableToPatchFieldRunTimeSelection\
(PatchTypeField, typePatchTypeField)                                           \
                                                                               \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// For non-templated patch fields
#define makePatchTypeField(PatchTypeField, typePatchTypeField)                 \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// For non-templated patch fields - use with caution
#define makeRemovablePatchTypeField(PatchTypeField, typePatchTypeField)        \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// For templated patch fields
#define makeTemplatePatchTypeField(fieldType, type)                            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );                                                                         \
    addToPatchFieldRunTimeSelection                                            \
    (                                                                          \
        CAT3(fvPatch, CAPITALIZE(fieldType), Field),                           \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field)                      \
    )

#define makePatchFields(type)                                                  \
    FOR_ALL_FIELD_TYPES(makeTemplatePatchTypeField, type)

#define makePatchFieldTypeName(fieldType, type)                                \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );

#define makePatchFieldTypeNames(type)                                          \
    FOR_ALL_FIELD_TYPES(makePatchFieldTypeName, type)

#define makePatchTypeFieldTypedef(fieldType, type)                             \
    typedef type##FvPatchField<fieldType>                                      \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field);

#define makePatchTypeFieldTypedefs(type)                                       \
    FOR_ALL_FIELD_TYPES(makePatchTypeFieldTypedef, type)


// error.H
//- Report an error message using Foam::FatalError
//  for functionName in file __FILE__ at line __LINE__
#define FatalErrorIn(functionName)                                             \
    ::Foam::FatalError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam::FatalError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)


//- Report an error message using Foam::FatalIOError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorIn(functionName, ios)                                      \
    ::Foam::FatalIOError((functionName), __FILE__, __LINE__, (ios))

//- Report an error message using Foam::FatalIOError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorInFunction(ios) FatalIOErrorIn(FUNCTION_NAME, ios)


//- Report an error message using Foam::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorIn(functionName, ios, msg)                             \
    ::Foam::IOerror::SafeFatalIOError                                          \
    ((functionName), __FILE__, __LINE__, (ios), (msg))

//- Report an error message using Foam::FatalIOError
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
        << "Not implemented" << ::Foam::abort(FatalError);

//- Issue a FatalErrorIn for a function not currently implemented.
//  The FUNCTION_NAME is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define NotImplemented notImplemented(FUNCTION_NAME)

// defineDebugSwitch.H
#define registerTemplateDebugSwitchWithName(Type,Name)                         \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<Type>                                      \
        Foam::RegisterDebugSwitch<Type>::registerDebugSwitch(Name)


//- Define the debug information, lookup as \a Name
#define registerDebugSwitchWithName(Type,Tag,Name)                             \
    class add##Tag##ToDebug                                                    \
    :                                                                          \
        public ::Foam::simpleRegIOobject                                       \
    {                                                                          \
    public:                                                                    \
        add##Tag##ToDebug(const char* name)                                    \
        :                                                                      \
            ::Foam::simpleRegIOobject(Foam::debug::addDebugObject, name)       \
        {}                                                                     \
        virtual ~add##Tag##ToDebug()                                           \
        {}                                                                     \
        virtual void readData(Foam::Istream& is)                               \
        {                                                                      \
            Type::debug = readLabel(is);                                       \
        }                                                                      \
        virtual void writeData(Foam::Ostream& os) const                        \
        {                                                                      \
            os << Type::debug;                                                 \
        }                                                                      \
    };                                                                         \
    add##Tag##ToDebug add##Tag##ToDebug_(Name)


//- Define the debug information, lookup as \a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch)                     \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch))

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch)                                   \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch);           \
    registerDebugSwitchWithName(Type, Type, Type::typeName_())

//- Define the debug information for templates, lookup as \a Name
#define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)             \
    template<>                                                                 \
    defineDebugSwitchWithName(Type, Name, DebugSwitch);                        \
    registerTemplateDebugSwitchWithName(Type, Name)

//- Define the debug information for templates sub-classes, lookup as \a Name
#define defineTemplate2DebugSwitchWithName(Type, Name, DebugSwitch)            \
    template<>                                                                 \
    defineDebugSwitchWithName(Type, Name, DebugSwitch);                        \
    registerTemplateDebugSwitchWithName(Type, Name)

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, DebugSwitch)                           \
    defineTemplateDebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, DebugSwitch)                      \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)


//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplate2DebugSwitch(Type, DebugSwitch)                          \
    defineTemplate2DebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplate2DebugSwitch(Type, DebugSwitch)                     \
    defineTemplate2DebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)


// className.H
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Declarations (without debug information)

//- Add typeName information from argument \a TypeNameString to a class.
//  Without debug information
#define ClassNameNoDebug(TypeNameString)                                       \
    static const char* typeName_() { return TypeNameString; }                  \
    static const ::Foam::word typeName

//- Add typeName information from argument \a TypeNameString to a namespace.
//  Without debug information.
#define NamespaceNameNoDebug(TypeNameString)                                   \
    inline const char* typeName_() { return TypeNameString; }                  \
    extern const ::Foam::word typeName

//- Add typeName information from argument \a TemplateNameString to a
//  template class.  Without debug information.
#define TemplateNameNoDebug(TemplateNameString)                                \
class TemplateNameString##Name                                                 \
{                                                                              \
public:                                                                        \
    TemplateNameString##Name() {}                                              \
    ClassNameNoDebug(#TemplateNameString);                                     \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Declarations (with debug information)

//- Add typeName information from argument \a TypeNameString to a class.
//  Also declares debug information.
#define ClassName(TypeNameString)                                              \
    ClassNameNoDebug(TypeNameString);                                          \
    static int debug

//- Add typeName information from argument \a TypeNameString to a namespace.
//  Also declares debug information.
#define NamespaceName(TypeNameString)                                          \
    NamespaceNameNoDebug(TypeNameString);                                      \
    extern int debug

//- Add typeName information from argument \a TypeNameString to a
//  template class.  Also declares debug information.
#define TemplateName(TemplateNameString)                                       \
class TemplateNameString##Name                                                 \
{                                                                              \
public:                                                                        \
    TemplateNameString##Name() {}                                              \
    ClassName(#TemplateNameString);                                            \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Definitions (without debug information)

//- Define the typeName, with alternative lookup as \a Name
#define defineTypeNameWithName(Type, Name)                                     \
    const ::Foam::word Type::typeName(Name)

//- Define the typeName
#define defineTypeName(Type)                                                   \
    defineTypeNameWithName(Type, Type::typeName_())

//- Define the typeName as \a Name for template classes
#define defineTemplateTypeNameWithName(Type, Name)                             \
    template<>                                                                 \
    defineTypeNameWithName(Type, Name)
//- Define the typeName as \a Name for template sub-classes
#define defineTemplate2TypeNameWithName(Type, Name)                            \
    template<>                                                                 \
    defineTypeNameWithName(Type, Name)

//- Define the typeName for template classes, useful with typedefs
#define defineTemplateTypeName(Type)                                           \
    defineTemplateTypeNameWithName(Type, #Type)

//- Define the typeName directly for template classes
#define defineNamedTemplateTypeName(Type)                                      \
    defineTemplateTypeNameWithName(Type, Type::typeName_())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Definitions (with debug information)

//- Define the typeName and debug information
#define defineTypeNameAndDebug(Type, DebugSwitch)                              \
    defineTypeName(Type);                                                      \
    defineDebugSwitch(Type, DebugSwitch)

//- Define the typeName and debug information, lookup as \a Name
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)        \
    defineTemplateTypeNameWithName(Type, Name);                                \
    defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)

//- Define the typeName and debug information for templates, useful
//  with typedefs
#define defineTemplateTypeNameAndDebug(Type, DebugSwitch)                      \
    defineTemplateTypeNameAndDebugWithName(Type, #Type, DebugSwitch)

//- Define the typeName and debug information for templates
#define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                 \
    defineNamedTemplateTypeName(Type);                                         \
    defineNamedTemplateDebugSwitch(Type, DebugSwitch)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// For templated sub-classes

//- Define the typeName and debug information, lookup as \a Name
#define defineTemplate2TypeNameAndDebugWithName(Type, Name, DebugSwitch)       \
    defineTemplate2TypeNameWithName(Type, Name);                               \
    defineTemplate2DebugSwitchWithName(Type, Name, DebugSwitch)

//- Define the typeName and debug information for templates, useful
//  with typedefs
#define defineTemplate2TypeNameAndDebug(Type, DebugSwitch)                     \
    defineTemplate2TypeNameAndDebugWithName(Type, #Type, DebugSwitch)

