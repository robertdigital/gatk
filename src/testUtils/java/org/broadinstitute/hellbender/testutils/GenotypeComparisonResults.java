package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.Genotype;

import java.util.*;
import java.util.stream.Collectors;

public class GenotypeComparisonResults {

    private Genotype actual;
    private Genotype expected;
    private Set<GenotypeAttributeEnum> mismatchedAttributes = new HashSet<>();
    private Set<String> mismatchedExtendedAttributes = new HashSet<>();

    public GenotypeComparisonResults(Genotype actual, Genotype expected){
        this.actual = actual;
        this.expected = expected;
    }

    public boolean isMatch(){
        return mismatchedAttributes.isEmpty();
    }

    public void addMismatchedAttribute(GenotypeAttributeEnum attributeName){
        mismatchedAttributes.add(attributeName);
    }

    public Set<GenotypeAttributeEnum> getMismatchedAttributes(){
        return mismatchedAttributes;
    }

    public void addMismatchedExtendedAttribute(String attribute){
        mismatchedExtendedAttributes.add(attribute);
    }

    public void addMismatchedExtendedAttributes(Collection<String> attributes){
        mismatchedExtendedAttributes.addAll(attributes);
    }

    public Set<String> getMismatchedExtendedAttributes(){
        return mismatchedExtendedAttributes;
    }

    public String getResultStringConcise(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("Genotype actual(" + actual.getSampleName() + ") matches expected(" + expected.getSampleName() + ")");
        }
        else{
            resultString.append("Genotype actual(" + actual.getSampleName() + ") does not match expected(" + expected.getSampleName() + ")");
            resultString.append(" with mismatches in the following attributes: " + mismatchedAttributes.stream().map(attr -> attr.getName()).collect(Collectors.joining(",")));
        }
        return resultString.toString();
    }

    public String getResultStringVerbose(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("Genotype actual(" + actual.getSampleName() + ") matches expected(" + expected.getSampleName() + ")");
        }
        else{
            resultString.append("Genotype actual(" + actual.getSampleName() + ") does not match expected(" + expected.getSampleName() + ") with mismatches in the following attributes:\n");
            Object actualValue = null;
            Object expectedValue = null;
            for(GenotypeAttributeEnum attribute : mismatchedAttributes){
                actualValue = attribute.getValue(actual);
                expectedValue = attribute.getValue(expected);
                resultString.append("  " + attribute.getName() + ": actual(" + actualValue + ") vs expected(" + expectedValue + ")\n");
                if(attribute.equals(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES) && !this.mismatchedExtendedAttributes.isEmpty()){
                    for(String extendedAttribute : this.mismatchedExtendedAttributes){
                        resultString.append("    " + extendedAttribute + ": actual(" + actual.getExtendedAttribute(extendedAttribute) + ") vs expected(" + expected.getExtendedAttribute(extendedAttribute) + ")\n");
                    }
                }
            }
        }
        return resultString.toString();
    }
}
