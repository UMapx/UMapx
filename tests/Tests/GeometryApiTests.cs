using System.ComponentModel;
using System.Drawing;
using System.Reflection;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class GeometryApiTests
{
    private const BindingFlags PublicMembers =
        BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly;

    private static Type Map(Type type)
    {
        if (type == typeof(Point)) return typeof(PointInt);
        if (type == typeof(PointF)) return typeof(PointFloat);
        if (type == typeof(Size)) return typeof(SizeInt);
        if (type == typeof(SizeF)) return typeof(SizeFloat);
        if (type == typeof(Rectangle)) return typeof(RectangleInt);
        if (type == typeof(RectangleF)) return typeof(RectangleFloat);
        if (type.IsArray) return Map(type.GetElementType()!).MakeArrayType();
        if (type.IsGenericType)
            return type.GetGenericTypeDefinition().MakeGenericType(type.GetGenericArguments().Select(Map).ToArray());
        return type;
    }

    private static void SameParameters(ParameterInfo[] expected, ParameterInfo[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.Equal(Map(expected[i].ParameterType), actual[i].ParameterType);
            Assert.Equal(expected[i].Name, actual[i].Name);
            Assert.Equal(expected[i].IsOptional, actual[i].IsOptional);
            Assert.Equal(expected[i].DefaultValue, actual[i].DefaultValue);
            Assert.Equal(expected[i].IsDefined(typeof(ParamArrayAttribute)),
                actual[i].IsDefined(typeof(ParamArrayAttribute)));
        }
    }

    [Theory]
    [InlineData(typeof(Point))]
    [InlineData(typeof(PointF))]
    [InlineData(typeof(Rectangle))]
    [InlineData(typeof(RectangleF))]
    [InlineData(typeof(Size))]
    [InlineData(typeof(SizeF))]
    public void AllSystemDrawingPublicMembersHaveMatchingSignatures(Type original)
    {
        var actual = Map(original);
        Assert.True(actual.IsValueType);
        foreach (var contract in original.GetInterfaces())
            Assert.Contains(Map(contract), actual.GetInterfaces());
        foreach (var constructor in original.GetConstructors())
        {
            var parameters = constructor.GetParameters();
            var match = actual.GetConstructor(parameters.Select(p => Map(p.ParameterType)).ToArray());
            Assert.NotNull(match);
            SameParameters(parameters, match.GetParameters());
        }
        foreach (var method in original.GetMethods(PublicMembers))
        {
            var parameters = method.GetParameters();
            var match = actual.GetMethods(PublicMembers).SingleOrDefault(m =>
                m.Name == method.Name && m.IsStatic == method.IsStatic &&
                m.ReturnType == Map(method.ReturnType) &&
                m.GetParameters().Select(p => p.ParameterType)
                    .SequenceEqual(parameters.Select(p => Map(p.ParameterType))));
            Assert.True(match != null, method.ToString());
            SameParameters(parameters, match!.GetParameters());
        }
        foreach (var property in original.GetProperties(PublicMembers))
        {
            var match = actual.GetProperty(property.Name);
            Assert.NotNull(match);
            Assert.Equal(Map(property.PropertyType), match.PropertyType);
            Assert.Equal(property.CanRead, match.CanRead);
            Assert.Equal(property.CanWrite, match.CanWrite);
            Assert.Equal(property.GetCustomAttribute<BrowsableAttribute>()?.Browsable,
                match.GetCustomAttribute<BrowsableAttribute>()?.Browsable);
        }
        foreach (var field in original.GetFields(PublicMembers))
        {
            var match = actual.GetField(field.Name);
            Assert.NotNull(match);
            Assert.Equal(Map(field.FieldType), match.FieldType);
            Assert.Equal(field.IsStatic, match.IsStatic);
            Assert.Equal(field.IsInitOnly, match.IsInitOnly);
            Assert.Equal(Activator.CreateInstance(actual), match.GetValue(null));
        }
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
    public void ConversionMethodsAndInstanceReadonlyModifiersAreAbsent(Type type)
    {
        foreach (var method in type.GetMethods(PublicMembers))
        {
            Assert.DoesNotContain(method.GetCustomAttributesData(),
                a => a.AttributeType.FullName == "System.Runtime.CompilerServices.IsReadOnlyAttribute");
            Assert.DoesNotContain(method.GetParameters(),
                p => p.ParameterType.Namespace == "System.Drawing");
            Assert.NotEqual("System.Drawing", method.ReturnType.Namespace);
        }
        if (type == typeof(RectangleInt) || type == typeof(RectangleFloat))
        {
            Assert.Null(type.GetMethod("op_Addition"));
            Assert.Null(type.GetMethod("op_Subtraction"));
        }
        if (type == typeof(PointFloat)) Assert.Null(type.GetMethod("Offset"));
    }
}
