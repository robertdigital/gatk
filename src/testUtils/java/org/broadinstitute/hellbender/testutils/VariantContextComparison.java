package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class VariantContextComparison {

    private VariantContext actual;
    private VariantContext expected;
    private Set<VariantContextAttributeEnum> attributesToIgnore;
    private List<String> extendedAttributesToIgnore;
    private Set<GenotypeAttributeEnum> genotypeAttributesToIgnore;
    private List<String> genotypeExtendedAttributesToIgnore;

    private VariantContextComparisonResults results;

    public static class Builder{
        private VariantContext actual;
        private VariantContext expected;
        private final Set<VariantContextAttributeEnum> variantContextAttributesToIgnore = new HashSet<>();
        private final List<String> variantContextExtendedAttributesToIgnore = new LinkedList<>();
        private final Set<GenotypeAttributeEnum> genotypeAttributesToIgnore = new HashSet<>();
        private final List<String> genotypeExtendedAttributesToIgnore = new LinkedList<>();

        public Builder(VariantContext actual, VariantContext expected){
            this.actual = actual;
            this.expected = expected;
        }

        public Builder addVariantContextAttributesToIgnore(Set<VariantContextAttributeEnum> variantContextAttributesToIgnore){
            this.variantContextAttributesToIgnore.addAll(variantContextAttributesToIgnore);

            return this;
        }

        public Builder addVariantContextAttributeToIgnore(VariantContextAttributeEnum variantContextAttributeToIgnore){
            this.variantContextAttributesToIgnore.add(variantContextAttributeToIgnore);

            return this;
        }

        public Builder addVariantContextExtendedAttributesToIgnore(List<String> variantContextExtendedAttributesToIgnore){
            this.variantContextExtendedAttributesToIgnore.addAll(variantContextExtendedAttributesToIgnore);

            return this;
        }

        public Builder addVariantContextExtendedAttributeToIgnore(String variantContextExtendedAttributeToIgnore){
            this.variantContextExtendedAttributesToIgnore.add(variantContextExtendedAttributeToIgnore);

            return this;
        }

        public Builder addGenotypeAttributesToIgnore(Set<GenotypeAttributeEnum> genotypeAttributesToIgnore){
            this.genotypeAttributesToIgnore.addAll(genotypeAttributesToIgnore);

            return this;
        }

        public Builder addGenotypeAttributeToIgnore(GenotypeAttributeEnum genotypeAttributeToIgnore){
            this.genotypeAttributesToIgnore.add(genotypeAttributeToIgnore);

            return this;
        }

        public Builder addGenotypeExtendedAttributesToIgnore(List<String> genotypeExtendedAttributesToIgnore){
            this.genotypeExtendedAttributesToIgnore.addAll(genotypeExtendedAttributesToIgnore);

            return this;
        }

        public Builder addGenotypeExtendedAttributeToIgnore(String genotypeExtendedAttributeToIgnore){
            this.genotypeExtendedAttributesToIgnore.add(genotypeExtendedAttributeToIgnore);

            return this;
        }

        public VariantContextComparison build(){
            VariantContextComparison comparison = new VariantContextComparison();
            comparison.actual = this.actual;
            comparison.expected = this.expected;
            comparison.attributesToIgnore = this.variantContextAttributesToIgnore;
            comparison.extendedAttributesToIgnore = this.variantContextExtendedAttributesToIgnore;
            comparison.genotypeAttributesToIgnore = this.genotypeAttributesToIgnore;
            comparison.genotypeExtendedAttributesToIgnore = this.genotypeExtendedAttributesToIgnore;

            return comparison;
        }
    }

    private VariantContextComparison(){};

    public VariantContext getActual() {
        return actual;
    }

    public void setActual(VariantContext actual) {
        this.actual = actual;
    }

    public VariantContext getExpected() {
        return expected;
    }

    public void setExpected(VariantContext expected) {
        this.expected = expected;
    }

    public Set<VariantContextAttributeEnum> getAttributesToIgnore() {
        return attributesToIgnore;
    }

    public void setAttributesToIgnore(Set<VariantContextAttributeEnum> attributesToIgnore) {
        this.attributesToIgnore = attributesToIgnore;
    }

    public List<String> getExtendedAttributesToIgnore() {
        return extendedAttributesToIgnore;
    }

    public void setExtendedAttributesToIgnore(List<String> extendedAttributesToIgnore) {
        this.extendedAttributesToIgnore = extendedAttributesToIgnore;
    }

    public Set<GenotypeAttributeEnum> getGenotypeAttributesToIgnore() {
        return genotypeAttributesToIgnore;
    }

    public void setGenotypeAttributesToIgnore(Set<GenotypeAttributeEnum> genotypeAttributesToIgnore) {
        this.genotypeAttributesToIgnore = genotypeAttributesToIgnore;
    }

    public List<String> getGenotypeExtendedAttributesToIgnore() {
        return genotypeExtendedAttributesToIgnore;
    }

    public void setGenotypeExtendedAttributesToIgnore(List<String> genotypeExtendedAttributesToIgnore) {
        this.genotypeExtendedAttributesToIgnore = genotypeExtendedAttributesToIgnore;
    }

    public VariantContextComparisonResults getResults(){
        return results;
    }

    public VariantContextComparison compare(){
        final VariantContextComparisonResults results = new VariantContextComparisonResults(this.actual, this.expected);

        if(actual == null || expected == null){
            throw new IllegalArgumentException("Variant contexts to be compared must not be null");
        }

        for(VariantContextAttributeEnum attribute : VariantContextAttributeEnum.values()){
            if(!attributesToIgnore.contains(attribute)){
                if(attribute == VariantContextAttributeEnum.ATTRIBUTES){
                    List<String> errorKeys = VariantContextTestUtils.checkAttributesEquals(
                            VariantContextTestUtils.filterIgnoredAttributes(actual.getAttributes(), extendedAttributesToIgnore),
                            VariantContextTestUtils.filterIgnoredAttributes(expected.getAttributes(), extendedAttributesToIgnore)
                    );
                    if(!errorKeys.isEmpty()){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.ATTRIBUTES);
                        results.addMismatchedExtendedAttributes(errorKeys);
                    }
                }
                else if(attribute == VariantContextAttributeEnum.GENOTYPES){
                    if((!actual.hasGenotypes() && expected.hasGenotypes()) || (actual.hasGenotypes() && !expected.hasGenotypes())){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                    }
                    else if(actual.hasGenotypes() && expected.hasGenotypes()){
                        for(String sampleName : expected.getSampleNames()){
                            if(actual.hasGenotype(sampleName)){
                                final GenotypeComparisonResults comparisonResults = new GenotypeComparison.Builder(actual.getGenotype(sampleName), expected.getGenotype(sampleName))
                                        .addAttributesToIgnore(this.genotypeAttributesToIgnore)
                                        .addExtendedAttributesToIgnore(this.genotypeExtendedAttributesToIgnore)
                                        .build()
                                        .compare().getResults();
                                results.addGenotypeComparisonResult(comparisonResults);
                                if(!results.getMismatchedAttributes().contains(VariantContextAttributeEnum.GENOTYPES) && !comparisonResults.isMatch()){
                                    results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                                }
                            }
                            else{
                                if(!results.getMismatchedAttributes().contains(VariantContextAttributeEnum.GENOTYPES)){
                                    results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                                }
                            }
                        }
                    }
                }
                else if(attribute == VariantContextAttributeEnum.PHRED_SCALED_QUALITY){
                    if(!BaseTest.equalsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual())){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.PHRED_SCALED_QUALITY);
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
