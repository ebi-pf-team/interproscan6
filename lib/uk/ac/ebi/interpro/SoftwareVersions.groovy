package uk.ac.ebi.interpro

import groovy.json.JsonOutput

class SoftwareVersions {
    static String toYaml(List entries) {
        def processes = new TreeMap<String, Map>()
        entries.each { entry ->
            def (process, tool, version) = entry
            def tools = processes.computeIfAbsent(process.toString()) { new TreeMap() }
            def versions = tools.computeIfAbsent(tool.toString()) { new TreeSet() }
            versions.add(version.toString().trim())
        }
        return processes.collect { process, tools ->
            def lines = tools.collect { tool, versions ->
                def value = versions.size() == 1 ? versions.first() : versions.toList()
                "  ${JsonOutput.toJson(tool)}: ${JsonOutput.toJson(value)}"
            }
            "${JsonOutput.toJson(process)}:\n${lines.join('\n')}\n"
        }.join('')
    }

    static String workflowVersion(String version, String commit) {
        return commit ? "${version}-${commit}" : version
    }
}
