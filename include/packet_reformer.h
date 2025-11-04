/**
 * @file packet_reformer.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class forreforming SCHP5 packets
 * @version 0.1
 * @date 2025-10-23
 */
# ifndef PACKET_REFORMER_H
# define PACKET_REFORMER_H

# include "simple_packet.h"
# include "crcinterface.h"
# include <cstdint>
# include <queue>

/** 
 * @brief Processes recieved SCHP5 protocol bitstream to packets
 */
class PacketReformer {

    private:

    /**
    * @brief The address of the sender
    */
    const uint8_t senderAddress;

    /**
    * @brief The address of the reciever
    */
    const uint8_t recieverAddress;

    /**
     * @brief If we are the data source (transaction role)
     */
    bool areWeSource;

    /**
     * @brief A pointer to the crc generator
     */
    crcutil_interface::CRC * crcGenny;

    /**
     * @brief A buffer to store control packets
     */
    std::vector<SimpPacket> * controlPacketBuffer;

    /**
     * @brief A buffer to store the current transaction's data packets
     */
    std::vector<SimpPacket> * dataPacketBuffer;
    
    /**
     * @brief The current processing position in the input data buffer
     */
    std::vector<bool>::iterator inputPos;

    public:

    /**
     * @brief The buffer of data to process into packets
     */
    std::vector<bool> * inDataBuffer;

    /**
     * @brief Construct a new packet reformer
     */
    PacketReformer(
        std::vector<bool> * in_data_buffer,
        std::vector<SimpPacket> * control_packet_buffer,
        std::vector<SimpPacket> * data_packet_buffer,
        uint8_t sender_address,
        uint8_t reciever_address,
        bool are_we_source);

    /**
     * @brief Tells if we are currently the source of data
     * @returns the current value of areWeSource
     */
    bool getAreWeSource();

    /**
     * @brief Sets whether we are curretly the source of data
     * @param are_we_source Whether we should be set as the source of data
     */
    void setAreWeSource(bool are_we_source);

    /**
     * @brief Gets the current processing position in the input data buffer
     * @returns The current position in the input data buffer
     */
    std::vector<bool>::iterator getInputPosition();

    /**
     * @brief Process raw bitstream to a packet
     * @returns 0 if not enough data in stream, 1 if bad due to missing end flag, 2 if for other reciever, 3 if from other sender, 254 if bad due to CRC, 255 if good
     */
    uint8_t formPacket();
    
    /**
     * @brief Processes packet in general buffer to either control or data packet buffers
     */
    void sortPacket(SimpPacket to_sort);
    
};

# endif