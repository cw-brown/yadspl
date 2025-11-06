/**
 * @file data_extractor.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for extracting the data from SCHP5 packets
 * @version 0.1
 * @date 2025-10-23
 */
# ifndef DATA_EXTRACTOR_H
# define DATA_EXTRACTOR_H

# include "simple_packet.h"
# include "crcinterface.h"
# include <cstdint>
# include <queue>

/** 
 * @brief Processes recieved SCHP5 protocol packets to raw bitstream
 */
class DataExtractor {

    private:

    /**
     * @brief A buffer to store which packets must be retransmitted
     */
    std::vector<uint8_t> * arqBuffer;

    /**
     * @brief A buffer to store which packets must be retransmitted
     */
    std::vector<uint8_t> * ackBuffer;

    /**
     * @brief A buffer to store the current transaction's data packets
     */
    std::vector<SimpPacket> * dataPacketBuffer;

    /**
     * @brief Store the packets that have been correctly recieved
     */
    std::vector<SimpPacket> * goodPackets;
    
    /**
     * @brief The current processing position in the input data buffer
     */
    std::vector<SimpPacket>::iterator inputPos;

    public:
    /**
     * @brief Construct a new data extractor
     */
    DataExtractor(
        std::vector<uint8_t> * arq_buffer,
        std::vector<uint8_t> * ack_buffer,
        std::vector<SimpPacket> * data_packet_buffer,
        std::vector<SimpPacket> * good_packets
        );

    /**
     * @brief Processes the next data packet in the buffer.
     */
    void processNextPacket();

    /**
     * @brief Gets the current processing position in the input data buffer
     * @returns The current position in the input data buffer
     */
    std::vector<SimpPacket>::iterator getInputPosition();
    
    /**
     * @brief Reports how many of the recieved packets are continuous
     * @returns The number of continuous good packets in the buffer, counting from the start
     */
    uint8_t getGoodPacketsContinuous();
    
    /**
     * @brief Reports whether the good packet buffer is full (all packets recieved in this transaction)
     * @returns True if the good packet buffer is full, false otherwise
     */
    bool getGoodPacketsFull();

    /**
     * @brief Sets all packets in the good packet buffer to the unrecieved state
     */
    void flushGoodPacketBuffer();
    
};

# endif