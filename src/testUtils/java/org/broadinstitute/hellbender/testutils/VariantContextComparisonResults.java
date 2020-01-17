package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;
import java.util.stream.Collectors;

public class VariantContextComparisonResults {

    private VariantContext actual;
    private VariantContext expected;
    private Set<VariantContextAttributeEnum> mismatchedAttributes = new HashSet<>();
    private Set<String> mismatchedExtendedAttributes = new HashSet<>();
    private List<GenotypeComparisonResults> genotypeComparisonResults = new LinkedList<>();


    public VariantContextComparisonResults(VariantContext actual, VariantContext expected){
        this.actual = actual;
        this.expected = expected;
    }

    public boolean isMatch(){
        return mismatchedAttributes.isEmpty();
    }

    public void addMismatchedAttribute(VariantContextAttributeEnum attributeName){
        mismatchedAttributes.add(attributeName);
    }

    public Set<VariantContextAttributeEnum> getMismatchedAttributes() {
        return mismatchedAttributes;
    }

    public void setMismatchedAttributes(Set<VariantContextAttributeEnum> mismatchedAttributes) {
        this.mismatchedAttributes = mismatchedAttributes;
    }

    public void addMismatchedExtendedAttribute(String attribute) {
        mismatchedExtendedAttributes.add(attribute);
    }

    public void addMismatchedExtendedAttributes(Collection<String> attributes) {
        mismatchedExtendedAttributes.addAll(attributes);
    }

    public Set<String> getMismatchedExtendedAttributes(){
        return mismatchedExtendedAttributes;
    }

    public void setMismatchedExtendedAttributes(Set<String> mismatchedExtendedAttributes) {
        this.mismatchedExtendedAttributes = mismatchedExtendedAttributes;
    }

    public void addGenotypeComparisonResult(GenotypeComparisonResults genotypeComparisonResult){
        this.genotypeComparisonResults.add(genotypeComparisonResult);
    }

    public List<GenotypeComparisonResults> getGenotypeComparisonResults() {
        return genotypeComparisonResults;
    }

    public void setGenotypeComparisonResults(List<GenotypeComparisonResults> genotypeComparisonResults) {
        this.genotypeComparisonResults = genotypeComparisonResults;
    }

    public String getResultStringConcise(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("VariantContext actual(" + actual.getSource() + ") matches expected(" + expected.getSource() + ")");
        }
        else{
            resultString.append("Variant actual(" + actual.getSource() + ") does not match expected(" + expected.getSource() + ")");
            resultString.append(" with mismatches in the following attributes: " + mismatchedAttributes.stream().map(attr -> attr.getName()).collect(Collectors.joining(",")));
        }
        return resultString.toString();
    }

    public String getResultStringVerbose(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("VariantContext actual(" + actual.getSource() + ") matches expected(" + expected.getSource() + ")");
        }
        else{
            resultString.append("VariantContext actual(" + actual.getSource() + ") does not match expected(" + expected.getSource() + ") with mismatches in the following attributes:\n");
            Object actualValue = null;
            Object expectedValue = null;
            for(VariantContextAttributeEnum attribute : mismatchedAttributes){
                actualValue = attribute.getValue(actual);
                expectedValue = attribute.getValue(expected);
                if(attribute.equals(VariantContextAttributeEnum.ATTRIBUTES) && !this.mismatchedExtendedAttributes.isEmpty()){
                    resultString.append("  " + attribute.getName() + ": actual(" + actualValue + ") vs expected(" + expectedValue + ")\n");
                    for(String extendedAttribute : this.mismatchedExtendedAttributes){
                        resultString.append("    " + extendedAttribute + ": actual(" + actual.getAttribute(extendedAttribute) + ") vs expected(" + expected.getAttribute(extendedAttribute) + ")\n");
                    }
                }
                else if(attribute.equals(VariantContextAttributeEnum.GENOTYPES)){
                    resultString.append("  Genotypes: actual(" + actual.getGenotypes().getSampleNamesOrderedByName().stream().collect(Collectors.joining(",")) + ") vs " +
                            "expected(" + expected.getGenotypes().getSampleNamesOrderedByName().stream().collect(Collectors.joining(",")) + ")\n");
                    for(GenotypeComparisonResults genotypeResults : genotypeComparisonResults){
                        genotypeResults.getResultStringVerbose();
                    }
                }
            }
        }
        return resultString.toString();
    }
}
