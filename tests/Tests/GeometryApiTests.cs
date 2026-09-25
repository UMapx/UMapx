using System.ComponentModel;
using System.Reflection;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public partial class GeometryApiTests
{
    private const BindingFlags PublicMembers =
        BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly;

    private static string TypeName(Type type) => type.IsGenericType
        ? type.Name.Split('`')[0] + "<" + string.Join(",", type.GetGenericArguments().Select(TypeName)) + ">"
        : type.Name;

    private static string Parameters(MethodBase method) => string.Join(", ", method.GetParameters().Select(p =>
        (p.IsDefined(typeof(ParamArrayAttribute)) ? "params " : "") + TypeName(p.ParameterType) + " " + p.Name
        + (p.IsOptional ? " optional" : "")));

    private static IEnumerable<string> Signatures(Type type)
    {
        foreach (var contract in type.GetInterfaces())
            yield return "interface " + TypeName(contract);
        foreach (var constructor in type.GetConstructors())
            yield return "new(" + Parameters(constructor) + ")";
        foreach (var method in type.GetMethods(PublicMembers).Where(m => !m.IsSpecialName || m.Name.StartsWith("op_")))
            yield return (method.IsStatic ? "static " : "") + TypeName(method.ReturnType) + " " + method.Name
                + "(" + Parameters(method) + ")";
        foreach (var property in type.GetProperties(PublicMembers))
            yield return "property " + ((property.GetMethod ?? property.SetMethod)!.IsStatic ? "static " : "")
                + TypeName(property.PropertyType) + " " + property.Name
                + (property.GetMethod?.IsPublic == true ? " get" : "")
                + (property.SetMethod?.IsPublic == true ? " set" : "")
                + (property.GetCustomAttribute<BrowsableAttribute>() is { } attribute ? " browsable=" + attribute.Browsable : "");
        foreach (var field in type.GetFields(PublicMembers))
            yield return "field " + (field.IsStatic ? "static " : "") + (field.IsInitOnly ? "readonly " : "")
                + TypeName(field.FieldType) + " " + field.Name;
    }

    [Theory]
    [MemberData(nameof(RequiredApi))]
    public void PublicMembersPreserveRequiredSignatures(Type type, string[] required)
    {
        Assert.True(type.IsValueType);
        Assert.Contains(typeof(ICloneable), type.GetInterfaces());
        var actual = Signatures(type).ToHashSet();
        foreach (var signature in required)
            Assert.Contains(signature, actual);
        Assert.Equal(Activator.CreateInstance(type), type.GetField("Empty")!.GetValue(null));
    }

    [Theory]
    [InlineData(typeof(PointInt))]
    [InlineData(typeof(PointFloat))]
    [InlineData(typeof(RectangleInt))]
    [InlineData(typeof(RectangleFloat))]
    [InlineData(typeof(SizeInt))]
    [InlineData(typeof(SizeFloat))]
    [InlineData(typeof(RangeInt))]
    [InlineData(typeof(RangeFloat))]
    public void ApiUsesNativeTypesWithoutInstanceReadonlyModifiers(Type type)
    {
        static void NativeType(Type value)
        {
            if (value.HasElementType) NativeType(value.GetElementType()!);
            else Assert.Contains(value.Namespace, new[] { "System", "System.Numerics", "UMapx.Core" });
        }
        foreach (var constructor in type.GetConstructors())
        foreach (var parameter in constructor.GetParameters())
            NativeType(parameter.ParameterType);
        foreach (var method in type.GetMethods(PublicMembers))
        {
            Assert.DoesNotContain(method.GetCustomAttributesData(),
                a => a.AttributeType.FullName == "System.Runtime.CompilerServices.IsReadOnlyAttribute");
            foreach (var parameter in method.GetParameters()) NativeType(parameter.ParameterType);
            NativeType(method.ReturnType);
        }
        if (type == typeof(RectangleInt) || type == typeof(RectangleFloat))
        {
            Assert.Null(type.GetMethod("op_Addition"));
            Assert.Null(type.GetMethod("op_Subtraction"));
        }
        if (type == typeof(PointFloat)) Assert.Null(type.GetMethod("Offset"));
    }
}
