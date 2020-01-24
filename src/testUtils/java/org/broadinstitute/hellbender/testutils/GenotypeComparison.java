package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.Genotype;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class GenotypeComparison {

    private Genotype actual;
    private Genotype expected;
    private Set<GenotypeAttributeEnum> attributesToIgnore;
    private List<String> extendedAttributesToIgnore;
    private boolean ignoreActualExtraAttributes;

    private GenotypeComparisonResults results;

    public static class Builder{

        private Genotype actual;
        private Genotype expected;
        private final Set<GenotypeAttributeEnum> attributesToIgnore = new HashSet<>();
        private final List<String> extendedAttributesToIgnore = new LinkedList<>();
        private boolean ignoreActualExtraAttributes = false;

        public Builder(Genotype actual, Genotype expected){
            if(actual == null || expected == null){
                throw new IllegalArgumentException("Genotypes being compared cannot be null");
            }
            this.actual = actual;
            this.expected = expected;
        }

        public Builder addAttributesToIgnore(Set<GenotypeAttributeEnum> attributesToIgnore){
            this.attributesToIgnore.addAll(attributesToIgnore);

            return this;
        }

        public Builder addAttributeToIgnore(GenotypeAttributeEnum attributeToIgnore){
            this.attributesToIgnore.add(attributeToIgnore);

            return this;
        }

        public Builder addExtendedAttributesToIgnore(List<String> extendedAttributesToIgnore){
            this.extendedAttributesToIgnore.addAll(extendedAttributesToIgnore);

            return this;
        }

        public Builder addExtendedAttributeToIgnore(String extendedAttributeToIgnore){
            this.extendedAttributesToIgnore.add(extendedAttributeToIgnore);

            return this;
        }

        public Builder ignoreActualExtraAttributes(){
            this.ignoreActualExtraAttributes = true;

            return this;
        }

        public GenotypeComparison build(){
            GenotypeComparison comparison = new GenotypeComparison();
            comparison.actual = this.actual;
            comparison.expected = this.expected;
            comparison.attributesToIgnore = this.attributesToIgnore;
            comparison.extendedAttributesToIgnore = this.extendedAttributesToIgnore;
            comparison.ignoreActualExtraAttributes = this.ignoreActualExtraAttributes;

            return comparison;
        }
    }

    public Genotype getActual() {
        return actual;
    }

    public void setActual(Genotype actual) {
        this.actual = actual;
    }

    public Genotype getExpected() {
        return expected;
    }

    public void setExpected(Genotype expected) {
        this.expected = expected;
    }

    public Set<GenotypeAttributeEnum> getAttributesToIgnore() {
        return attributesToIgnore;
    }

    public void setAttributesToIgnore(Set<GenotypeAttributeEnum> attributesToIgnore) {
        this.attributesToIgnore = attributesToIgnore;
    }

    public List<String> getExtendedAttributesToIgnore() {
        return extendedAttributesToIgnore;
    }

    public void setExtendedAttributesToIgnore(List<String> extendedAttributesToIgnore) {
        this.extendedAttributesToIgnore = extendedAttributesToIgnore;
    }

    public GenotypeComparisonResults getResults(){
        return results;
    }

    public GenotypeComparison compare(){
        GenotypeComparisonResults results = new GenotypeComparisonResults(actual, expected);

        for(GenotypeAttributeEnum attribute : GenotypeAttributeEnum.values()){
            if(!attributesToIgnore.contains(attribute)){
                if(attribute == GenotypeAttributeEnum.EXTENDED_ATTRIBUTES){
                    List<String> errorKeys = VariantContextTestUtils.checkAttributesEquals(
                            VariantContextTestUtils.filterIgnoredAttributes(actual.getExtendedAttributes(), extendedAttributesToIgnore),
                            VariantContextTestUtils.filterIgnoredAttributes(expected.getExtendedAttributes(), extendedAttributesToIgnore),
                            ignoreActualExtraAttributes
                    );
                    if(!errorKeys.isEmpty()){
                        results.addMismatchedAttribute(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES);
                        results.addMismatchedExtendedAttributes(errorKeys);
                    }
                }
                else{
                    if(!VariantContextTestUtils.checkFieldEquals(attribute.getValue(actual), attribute.getValue(expected))){
                        results.addMismatchedAttribute(attribute);
                    }
                }
            }
        }

        this.results = results;
        return this;
    }
}
