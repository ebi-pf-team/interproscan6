package uk.ac.ebi.interpro

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption

class LocalSeqDB {
    static Path create(String fileName, Object configuredDir) {
        Path localDir = createTempDir(configuredDir)
        return localDir.resolve(fileName)
    }

    static Path copyFrom(Path source, Object configuredDir) {
        Path localPath = create(source.fileName.toString(), configuredDir)
        Files.copy(source, localPath, StandardCopyOption.REPLACE_EXISTING)
        return localPath
    }

    static void copyTo(Path source, Path target) {
        Files.copy(source, target, StandardCopyOption.REPLACE_EXISTING)
    }

    static void cleanup(Path dbPath) {
        if (dbPath == null) {
            return
        }
        ["", "-wal", "-shm", "-journal"].each { suffix ->
            Files.deleteIfExists(Paths.get(dbPath.toString() + suffix))
        }
        Files.deleteIfExists(dbPath.parent)
    }

    private static Path createTempDir(Object configuredDir) {
        Path baseDir = getBaseDir(configuredDir)
        Files.createDirectories(baseDir)
        return Files.createTempDirectory(baseDir, "ips6-seqdb-")
    }

    private static Path getBaseDir(Object configuredDir) {
        if (configuredDir != null && configuredDir.toString().trim()) {
            return Paths.get(configuredDir.toString())
        }
        String envTmp = System.getenv("TMPDIR") ?: System.getenv("TMP") ?: System.getenv("TEMP")
        return Paths.get(envTmp ?: System.getProperty("java.io.tmpdir"))
    }
}
