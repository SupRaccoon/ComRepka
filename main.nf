#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { COMPARISON } from './workflows/comparison'

workflow {
    COMPARISON()
}
