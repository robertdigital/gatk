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

    private GenotypeComparisonResults results;

    public static class Builder{

        private Genotype actual;
        private Genotype expected;
        private final Set<GenotypeAttributeEnum> attributesToIgnore = new HashSet<>();
        private final List<String> extendedAttributesToIgnore = new LinkedList<>();

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

        public GenotypeComparison build(){
            GenotypeComparison comparison = new GenotypeComparison();
            comparison.actual = this.actual;
            comparison.expected = this.expected;
            comparison.attributesToIgnore = this.attributesToIgnore;
            comparison.extendedAttributesToIgnore = this.extendedAttributesToIgnore;

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
                            VariantContextTestUtils.filterIgnoredAttributes(expected.getExtendedAttributes(), extendedAttributesToIgnore)
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
        /* TODO: Remove this
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.SAMPLE_NAME)){
            if(!VariantContextTestUtils.checkFieldEquals(actual.getSampleName(), expected.getSampleName())){
                results.addMismatchedAttribute(GenotypeAttributeEnum.SAMPLE_NAME);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.ALLELES)){
            if(!VariantContextTestUtils.checkFieldEquals(actual.getAlleles(), expected.getAlleles())){
                results.addMismatchedAttribute(GenotypeAttributeEnum.ALLELES);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.GENOTYPE_STRING)){
            if(!VariantContextTestUtils.checkFieldEquals(actual.getGenotypeString(), expected.getGenotypeString())){
                results.addMismatchedAttribute(GenotypeAttributeEnum.GENOTYPE_STRING);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.TYPE)){
            if(!VariantContextTestUtils.checkFieldEquals(actual.getType(), expected.getType())){
                results.addMismatchedAttribute(GenotypeAttributeEnum.TYPE);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.FILTERS)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getFilters(), expected.getFilters())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.FILTERS);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.DP)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getDP(), expected.getDP())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.DP);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.AD)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getAD(), expected.getAD())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.AD);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.GQ)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getGQ(), expected.getGQ())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.GQ);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.PL)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getPL(), expected.getPL())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.PL);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.LIKELIHOODS)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getLikelihoods(), expected.getLikelihoods())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.LIKELIHOODS);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.IS_PHASED)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.isPhased(), expected.isPhased())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.IS_PHASED);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.PLOIDY)) {
            if (!VariantContextTestUtils.checkFieldEquals(actual.getPloidy(), expected.getPloidy())) {
                results.addMismatchedAttribute(GenotypeAttributeEnum.PLOIDY);
            }
        }
        if(!attributesToIgnore.contains(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES)) {
            List<String> errorKeys = VariantContextTestUtils.checkAttributesEquals(
                    VariantContextTestUtils.filterIgnoredAttributes(actual.getExtendedAttributes(), extendedAttributesToIgnore),
                    VariantContextTestUtils.filterIgnoredAttributes(expected.getExtendedAttributes(), extendedAttributesToIgnore)
            );
            if(!errorKeys.isEmpty()){
                results.addMismatchedAttribute(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES);
                results.addMismatchedExtendedAttributes(errorKeys);
            }
        }
        */

        this.results = results;
        return this;
    }
}
