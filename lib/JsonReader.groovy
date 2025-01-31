import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonParseException
import com.fasterxml.jackson.databind.JsonMappingException
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper

class JsonReader {
    static stream(String filePath, ObjectMapper mapper, Closure closure) {
        /*
        Stream the data, one jsonNode at a time from a Json file containing a JsonArray.
        It requires less memory to store maps as JsonNodes than as Maps,
        especially when the maps contain nested maps.
        To use:
        ObjectMapper mapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
        JsonReader.stream(filePath, mapper) { JsonNode node ->
            md5 = node.get("md5").asText()
        }
         */
        try {
            JsonFactory factory = mapper.getFactory()
            JsonParser parser = factory.createParser(new File(filePath))
            if (parser.nextToken() == JsonToken.START_ARRAY) {
                while (parser.nextToken() != JsonToken.END_ARRAY) {
                    JsonNode node = mapper.readTree(parser)
                    closure.call(node)
                }
            }
            parser.close()
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $filePath", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath", e)
        }
    }

    static Map<String, JsonNode> load(String filePath, ObjectMapper mapper) {
        /*
        Load all contents of a json file into the memory, storing internal nodes as JsonNodes
        to reduce memory requirements compared to storing the data in memory as maps.
        To use:
        ObjectMapper mapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
        def matchesMap = JsonReader.load(filePath, mapper)
        seqMatches = matchesMap.get("seqId_1")
        */
        try {
            File file = new File(filePath)
            JsonNode rootNode = mapper.readTree(file)
            if (!rootNode.isObject()) {
                throw new JsonMappingException("Expected JSON object at root level")
            }

            Map<String, JsonNode> dataMap = new HashMap<>()
            rootNode.fields().each { entry ->
                dataMap.put(entry.key, entry.value)
            }

            return dataMap
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $filePath", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath", e)
        }
    }
}
